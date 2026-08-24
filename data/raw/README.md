# Raw source data (not tracked in Git)

These files exceed GitHub's file-size limits and are the immutable Kaggle source, so they
are not stored here. Obtain them one of two ways:

1. Kaggle — "Structural Protein Sequences"
   https://www.kaggle.com/datasets/shahir/protein-data-set
   Snapshot used: [AUTHOR: record the snapshot date/version]
   Files: pdb_data_seq.csv (467,304 chains: structureId, chainId, sequence, residueCount,
   classification) and pdb_data_no_dups.csv (141,401 structure-level metadata rows).

2. Zenodo — the archived snapshot deposited with this project (definitive input):
   https://doi.org/[AUTHOR: insert Zenodo DOI]

Place both CSVs in this directory before running the pipeline from raw. To skip this step,
use data/processed/dataset_clustered.csv, which is included in the repository.
