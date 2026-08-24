# Intermediate outputs (not tracked in Git)

Regenerable search chunks produced during the run and archived on Zenodo rather than in Git:

- uniref50_chunks_fast/ , uniref50_chunks_sens/ — tiered DIAMOND UniRef50 search chunks (Section 6d)
- pretraining_strat_chunks/ — per-seed two-way stratification outputs (Section 6e)
- mmseqs/ — MMseqs2 scratch (matrices, temporary databases)

All of these are recreated deterministically by the pipeline notebook; they are archived
only to allow exact re-inspection without re-running the searches.
