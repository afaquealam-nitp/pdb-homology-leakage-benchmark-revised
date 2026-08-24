# Data

## processed/dataset_clustered.csv  (included in Git, ~13 MB)
The curated benchmark: 40,211 unique protein sequences across ten functional classes,
with 30%-identity cluster assignments. Columns: `structureId, chainId, sequence,
classification, L, id, cluster`. This is the frozen input to every experiment — the
identity-controlled evaluation reproduces from this file alone.

## raw/  (NOT in Git — download or fetch from Zenodo)
Source data redistributed from the RCSB Protein Data Bank via Kaggle. See raw/README.md.
`dataset_clustered.csv` is produced from these by the pipeline notebook.

## Provenance chain
Kaggle raw (pdb_data_seq.csv + pdb_data_no_dups.csv)
  → label audit, dedup, length/canonical filters, multi-label removal, top-10 classes
  → 30% MMseqs2 clustering
  → data/processed/dataset_clustered.csv (40,211 sequences)

## Checksums
Record SHA-256 for each raw file here before archiving, e.g.:
    sha256sum data/raw/*.csv > data/CHECKSUMS.txt
