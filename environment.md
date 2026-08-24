# Environment and tool versions

Fill in the exact versions used before archiving (this file is the record referenced
by the manuscript's Data Availability statement).

## Python packages
Run `pip freeze > requirements-lock.txt` in the environment used for the final run and
commit it. Key packages: numpy, pandas, scikit-learn, scipy, torch, transformers,
matplotlib, tqdm.

## Protein language model
- Checkpoint: `esm2_t33_650M_UR50D` (loaded via Hugging Face `transformers`)

## External command-line tools
- MMseqs2 — clustering (`easy-cluster`) and alignment top-hit transfer. Version: [AUTHOR: fill in]
- DIAMOND — UniRef50 exposure search (`blastp`, tiered sensitivity). Version: [AUTHOR: fill in]

## Evaluation
- Seeds: 0–9 (ten seeds)
- Hardware: single GPU for ESM-2 embedding; embeddings cached. [AUTHOR: note GPU model]
