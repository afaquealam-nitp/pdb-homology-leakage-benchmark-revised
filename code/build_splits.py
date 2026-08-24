#!/usr/bin/env python3
"""
build_splits.py — Homology-controlled train/val/test splits + leakage audit.

Companion to eval/metrics.py for the HySTraG-IP contact-refinement benchmark.

What it does
------------
1. Clusters all chain sequences with MMseqs2 at a chosen identity (default 30%).
2. Assigns *whole clusters* to train/val/test (a cluster is never split), which is
   the key step that prevents near-duplicate leakage.
3. Runs a leakage audit: for every test chain, the maximum sequence identity to any
   training chain, and the fraction of test chains with any hit above the clustering
   threshold. A clean split drives this near zero.
4. (Optional) Produces a temporal split by PDB deposition date as a robustness check.

Outputs (in --out_dir)
----------------------
  train.txt / val.txt / test.txt   one chain id per line
  splits.json                      {chain_id: split}
  clusters.tsv                     representative<TAB>member
  leakage_audit.csv                test_chain, max_identity_to_train, hit_above_threshold
  leakage_summary.md               a small table to paste into the manuscript

Requirements
------------
  * MMseqs2 on PATH (https://github.com/soedinglab/MMseqs2). No Python deps beyond stdlib.
  * Input FASTA whose headers are the chain ids used everywhere else in the pipeline,
    e.g.  >1ABC_A
  * (temporal mode) a CSV with columns: chain_id, deposition_date  (YYYY-MM-DD)

Example
-------
  python build_splits.py --fasta all_chains.fasta --out_dir splits/ \
      --min-seq-id 0.30 --coverage 0.8 --split 0.8 0.1 0.1 --seed 0
  python build_splits.py --fasta all_chains.fasta --out_dir splits_time/ \
      --mode temporal --metadata chains.csv --cutoff-date 2021-05-01
"""
import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


# --------------------------------------------------------------------------- #
# FASTA helpers
# --------------------------------------------------------------------------- #
def read_fasta(path):
    """Return dict {chain_id: sequence}. Header id = first whitespace token, '>' stripped."""
    seqs, cid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cid is not None:
                    seqs[cid] = "".join(buf)
                cid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cid is not None:
        seqs[cid] = "".join(buf)
    return seqs


def write_fasta(seqs, ids, path):
    with open(path, "w") as fh:
        for cid in ids:
            fh.write(f">{cid}\n{seqs[cid]}\n")


# --------------------------------------------------------------------------- #
# MMseqs2 wrappers
# --------------------------------------------------------------------------- #
def _require(tool):
    if shutil.which(tool) is None:
        sys.exit(f"[error] '{tool}' not found on PATH. Install it and try again.")


def mmseqs_cluster(fasta, out_prefix, tmp, min_seq_id, coverage, mmseqs="mmseqs"):
    """Cluster with `mmseqs easy-cluster`. Returns dict {member: representative}."""
    _require(mmseqs)
    cmd = [
        mmseqs, "easy-cluster", fasta, out_prefix, tmp,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", "1",
    ]
    subprocess.run(cmd, check=True)
    member2rep = {}
    with open(out_prefix + "_cluster.tsv") as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            member2rep[member] = rep
    return member2rep


def mmseqs_max_identity(query_fasta, target_fasta, out_dir, tmp, mmseqs="mmseqs"):
    """For each query chain, the max fraction-identity to any target chain.
    Uses `mmseqs easy-search`. Returns dict {query_id: max_fident} (0.0 if no hit)."""
    _require(mmseqs)
    aln = os.path.join(out_dir, "_leak_aln.m8")
    cmd = [
        mmseqs, "easy-search", query_fasta, target_fasta, aln, tmp,
        "-s", "7.5", "--max-seqs", "300",
        "--format-output", "query,target,fident,alnlen,qcov,tcov",
    ]
    subprocess.run(cmd, check=True)
    best = defaultdict(float)
    with open(aln) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, _t, fident = parts[0], parts[1], float(parts[2])
            if q == _t:          # skip self-hit if query happens to also be in target
                continue
            if fident > best[q]:
                best[q] = fident
    return best


# --------------------------------------------------------------------------- #
# Splitting
# --------------------------------------------------------------------------- #
def cluster_split(member2rep, fractions, seed):
    """Assign whole clusters to splits. Largest-cluster-first greedy toward target counts.
    fractions: (train, val, test) summing to ~1. Returns dict {chain_id: split}."""
    clusters = defaultdict(list)
    for member, rep in member2rep.items():
        clusters[rep].append(member)

    total = sum(len(m) for m in clusters.values())
    targets = {"train": fractions[0] * total,
               "val": fractions[1] * total,
               "test": fractions[2] * total}
    counts = {"train": 0, "val": 0, "test": 0}
    assign = {}

    # deterministic shuffle, then sort by size desc so big clusters are placed first
    reps = list(clusters.keys())
    random.Random(seed).shuffle(reps)
    reps.sort(key=lambda r: len(clusters[r]), reverse=True)

    for rep in reps:
        members = clusters[rep]
        # pick the split furthest below its target (by remaining headroom)
        split = max(("train", "val", "test"), key=lambda s: targets[s] - counts[s])
        for m in members:
            assign[m] = split
        counts[split] += len(members)
    return assign, counts, len(clusters)


def temporal_split(seqs, metadata_csv, cutoff_date, val_frac, seed):
    """train = deposited < cutoff; {val,test} = deposited >= cutoff (split by chain, seeded).
    metadata_csv columns: chain_id, deposition_date (YYYY-MM-DD)."""
    date = {}
    with open(metadata_csv) as fh:
        for row in csv.DictReader(fh):
            date[row["chain_id"]] = row["deposition_date"]
    missing = [c for c in seqs if c not in date]
    if missing:
        print(f"[warn] {len(missing)} chains lack a deposition_date; excluded from temporal split.")
    train = [c for c in seqs if date.get(c, "9999") < cutoff_date]
    after = [c for c in seqs if date.get(c, "0000") >= cutoff_date and c in date]
    random.Random(seed).shuffle(after)
    n_val = int(round(val_frac * len(after)))
    val, test = after[:n_val], after[n_val:]
    assign = {}
    for c in train: assign[c] = "train"
    for c in val:   assign[c] = "val"
    for c in test:  assign[c] = "test"
    return assign, {"train": len(train), "val": len(val), "test": len(test)}, None


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def write_outputs(out_dir, seqs, assign, member2rep):
    os.makedirs(out_dir, exist_ok=True)
    for split in ("train", "val", "test"):
        ids = sorted(c for c, s in assign.items() if s == split)
        with open(os.path.join(out_dir, f"{split}.txt"), "w") as fh:
            fh.write("\n".join(ids) + ("\n" if ids else ""))
    with open(os.path.join(out_dir, "splits.json"), "w") as fh:
        json.dump(assign, fh, indent=0, sort_keys=True)
    if member2rep:
        with open(os.path.join(out_dir, "clusters.tsv"), "w") as fh:
            for member, rep in sorted(member2rep.items()):
                fh.write(f"{rep}\t{member}\n")


def leakage_report(out_dir, seqs, assign, threshold, tmp, mmseqs):
    train_ids = [c for c, s in assign.items() if s == "train"]
    test_ids = [c for c, s in assign.items() if s == "test"]
    if not train_ids or not test_ids:
        print("[warn] leakage audit skipped (empty train or test).")
        return None
    tr_fa = os.path.join(out_dir, "_train.fasta")
    te_fa = os.path.join(out_dir, "_test.fasta")
    write_fasta(seqs, train_ids, tr_fa)
    write_fasta(seqs, test_ids, te_fa)
    best = mmseqs_max_identity(te_fa, tr_fa, out_dir, tmp, mmseqs)

    rows = [(c, round(best.get(c, 0.0), 4), int(best.get(c, 0.0) > threshold)) for c in sorted(test_ids)]
    with open(os.path.join(out_dir, "leakage_audit.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["test_chain", "max_identity_to_train", "hit_above_threshold"])
        w.writerows(rows)

    n = len(rows)
    n_hit = sum(r[2] for r in rows)
    ids = [r[1] for r in rows]
    ids_sorted = sorted(ids)
    median = ids_sorted[n // 2] if n else 0.0
    mx = max(ids) if ids else 0.0
    pct = 100.0 * n_hit / n if n else 0.0
    with open(os.path.join(out_dir, "leakage_summary.md"), "w") as fh:
        fh.write("| Leakage audit (test vs. train) | Value |\n|---|---|\n")
        fh.write(f"| Test chains | {n} |\n")
        fh.write(f"| Median max identity to train | {median:.3f} |\n")
        fh.write(f"| Maximum max identity to train | {mx:.3f} |\n")
        fh.write(f"| Test chains with a hit > {threshold:.2f} | {n_hit} ({pct:.1f}%) |\n")
    return dict(n=n, n_hit=n_hit, pct=pct, median=median, max=mx)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="FASTA of all chains; headers are chain ids.")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--mode", choices=["cluster", "temporal"], default="cluster")
    ap.add_argument("--min-seq-id", type=float, default=0.30, help="MMseqs2 clustering identity.")
    ap.add_argument("--coverage", type=float, default=0.8)
    ap.add_argument("--split", type=float, nargs=3, default=[0.8, 0.1, 0.1],
                    metavar=("TRAIN", "VAL", "TEST"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--metadata", help="CSV: chain_id,deposition_date (temporal mode).")
    ap.add_argument("--cutoff-date", help="YYYY-MM-DD (temporal mode).")
    ap.add_argument("--val-frac", type=float, default=0.5,
                    help="Fraction of post-cutoff chains used for val (temporal mode).")
    ap.add_argument("--leakage-threshold", type=float, default=None,
                    help="Identity above which a test->train hit counts as leakage "
                         "(default = --min-seq-id).")
    ap.add_argument("--no-audit", action="store_true", help="Skip the leakage audit.")
    ap.add_argument("--mmseqs", default="mmseqs")
    args = ap.parse_args()

    if abs(sum(args.split) - 1.0) > 1e-6:
        sys.exit("[error] --split must sum to 1.0")
    os.makedirs(args.out_dir, exist_ok=True)
    # MMseqs2 executes a shell script in its temp dir; a Google Drive FUSE mount forbids
    # executing files (error 13). Keep tmp on LOCAL disk even if out_dir is on Drive.
    tmp = os.path.join(tempfile.gettempdir(), "mmseqs_tmp")
    thr = args.leakage_threshold if args.leakage_threshold is not None else args.min_seq_id

    seqs = read_fasta(args.fasta)
    print(f"[info] {len(seqs)} chains read from {args.fasta}")

    member2rep = {}
    if args.mode == "cluster":
        prefix = os.path.join(args.out_dir, "clust")
        member2rep = mmseqs_cluster(args.fasta, prefix, tmp, args.min_seq_id,
                                    args.coverage, args.mmseqs)
        assign, counts, n_clust = cluster_split(member2rep, args.split, args.seed)
        print(f"[info] {n_clust} clusters -> "
              f"train {counts['train']} / val {counts['val']} / test {counts['test']}")
    else:
        if not (args.metadata and args.cutoff_date):
            sys.exit("[error] temporal mode needs --metadata and --cutoff-date")
        assign, counts, _ = temporal_split(seqs, args.metadata, args.cutoff_date,
                                           args.val_frac, args.seed)
        print(f"[info] temporal split @ {args.cutoff_date} -> "
              f"train {counts['train']} / val {counts['val']} / test {counts['test']}")

    write_outputs(args.out_dir, seqs, assign, member2rep)

    if not args.no_audit:
        summary = leakage_report(args.out_dir, seqs, assign, thr, tmp, args.mmseqs)
        if summary:
            print(f"[info] leakage audit: {summary['n_hit']}/{summary['n']} test chains "
                  f"({summary['pct']:.1f}%) have a >|{thr:.2f}| identity hit to train; "
                  f"median={summary['median']:.3f}, max={summary['max']:.3f}")
    print(f"[done] wrote splits and reports to {args.out_dir}")


if __name__ == "__main__":
    main()
