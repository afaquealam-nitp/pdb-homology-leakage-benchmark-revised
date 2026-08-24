# Reproducing the results

## Option A — from the curated dataset (fast, no raw download)
1. `pip install -r requirements.txt` and install MMseqs2 + DIAMOND (see ../environment.md).
2. Open `code/pdb_leakage_rebenchmark_sweep_6models.ipynb`.
3. Point it at `data/processed/dataset_clustered.csv`.
4. Run the identity-controlled evaluation (random vs cluster split, ten seeds) and the
   threshold sweep. Outputs are written to `results/tables/` and `results/figures/`.

## Option B — full pipeline from raw
1. Download the two raw CSVs into `data/raw/` (see data/raw/README.md).
2. Run the dataset-construction cells (label audit → dedup → filters → top-10 classes →
   30% MMseqs2 clustering) to regenerate `data/processed/dataset_clustered.csv`.
3. Continue as in Option A.

## Pretraining-exposure analysis (Sections 6d–6e)
Requires the UniRef50 database for the DIAMOND search. The tiered-sensitivity search and the
two-way (training-identity × pretraining-proxy-identity) stratification write
`uniref50_exposure.csv`, `pretraining_stratified_results.csv`, and
`pretraining_stratified_summary.csv`.

## Notes
- All randomness is seeded (0–9); report mean ± s.d. / 95% CI over seeds.
- ESM-2 embeddings are computed once and cached; the frozen probe/1-NN are deterministic given the cache.
- Runtime is dominated by the ESM-2 embedding pass and the UniRef50 search.
