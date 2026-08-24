# PDB Homology-Leakage Benchmark

Code, data, and results for the paper:

Afaque Alam and Mukesh Kumar, Department of Computer Science and Engineering, National Institute of Technology Patna.


## What this is

An audit of homology leakage on a widely used protein functional-class benchmark — the RCSB "Structural Protein Sequences" collection redistributed through Kaggle. From 467,304 chains we build a label-audited set of **40,211 unique sequences** across **ten functional classes** (**9,533 clusters** at 30% identity) and evaluate six models under a conventional random split and under identity-controlled cluster splits at MMseqs2 thresholds from 90% down to 20%, with ten seeds and confidence intervals.

Six models: a 1-NN k-mer memorization probe, Linear SVM, Logistic Regression, MMseqs2 alignment top-hit transfer, a frozen ESM-2 650M 1-NN, and a frozen ESM-2 650M linear probe (`esm2_t33_650M_UR50D`).

Headline findings: 86.5% of random-split test sequences have a training homolog above 30% identity; composition-based models lose roughly two-thirds of their apparent performance once homology is controlled; and the ranking of alignment transfer vs. the language-model representation inverts near 40% identity. A two-way stratification (identity to the supervised training set × identity to a UniRef50 pretraining proxy) shows the doubly-distant stratum needed to separate generalization from pretraining exposure is essentially empty, so this benchmark cannot settle that question.

## Repository layout

```
code/        pipeline notebook + build_splits.py + metrics.py
data/
  raw/       Kaggle source (NOT in Git — see data/raw/README.md)
  processed/ dataset_clustered.csv — the 40,211-sequence curated benchmark
results/
  tables/    all result CSVs (per-seed scores, sweep, exposure, stratification)
  figures/   fig1–fig5 as in the manuscript
  interim/   regenerable search chunks (NOT in Git — archived on Zenodo)
docs/        data-availability statement + reproduction guide
```

## Quick start

```bash
pip install -r requirements.txt
# install external tools (see environment.md): MMseqs2, DIAMOND
jupyter notebook code/pdb_leakage_rebenchmark_sweep_6models.ipynb
```

`data/processed/dataset_clustered.csv` is the frozen input to every experiment and is included here, so the identity-controlled evaluation reproduces without re-downloading the raw Kaggle files. See `docs/REPRODUCE.md` for the full pipeline.

## Data and archival

The large raw inputs and regenerable intermediate chunks are not stored in Git (GitHub file-size limits). They are archived on Zenodo alongside a complete snapshot of this repository; see `docs/DATA_AVAILABILITY.md`. Raw source: Kaggle "Structural Protein Sequences" (https://www.kaggle.com/datasets/shahir/protein-data-set), derived from the RCSB Protein Data Bank.

## Citation

See `CITATION.cff`. Code is released under MIT (`LICENSE`); data and results under CC-BY-4.0 (`LICENSE-DATA`).
