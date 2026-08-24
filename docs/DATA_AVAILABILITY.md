# Data availability

The source data are publicly available from Kaggle
(https://www.kaggle.com/datasets/shahir/protein-data-set), derived from the RCSB Protein
Data Bank. The curated dataset, cluster assignments, per-seed results, and the complete
reproducible pipeline are in this repository
(https://github.com/afaquealam-nitp/pdb-homology-leakage-benchmark).

A permanently archived snapshot — code, curated dataset, results, figures, and the raw and
intermediate data — is deposited at Zenodo: https://doi.org/[AUTHOR: insert Zenodo DOI].

For exact reproduction we pin: the Kaggle snapshot ([AUTHOR: date/version]); the protein
language-model checkpoint (esm2_t33_650M_UR50D); the MMseqs2 and DIAMOND versions and all
clustering and search parameters ([AUTHOR: versions], also in environment.md); and the ten
evaluation seeds (0–9). Because Kaggle collections can change over time, the deposited Zenodo
snapshot rather than the live Kaggle link is the definitive input.
