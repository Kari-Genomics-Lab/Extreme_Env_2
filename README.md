# ExtPrime: Extremophile Genomic-Signature Experiments

**TL;DR:** Build composite genome proxies and run supervised learning experiments to classify extremophile organisms by temperature and pH tolerance using k-mer features.

This repository provides pipelines to build **composite genome proxies**, compute **k-mer features**, and run **supervised experiments** on extremophile datasets (Temperature & pH). Inputs live in `data/`; all generated results go to `outputs/`.

---

## Installation

```bash
# From repo root
pip install -e .
```
## Download Data
Download the extremophile genome assemblies and metadata from Zenodo DOI: [10.5281/zenodo.17148766]

## Data Layout

Place downloaded assemblies (FASTA .fna) and the metadat file under `data/`
All results will be written to `outputs/` automatically.
All the results of the experiments of this study are available in `results/` folder.

## Supervised Learning Experiments

### 1) Effect of Genome Proxy Selection (Multiple Runs)

Tests whether random genome proxy choice changes classification accuracy.

```bash
python3 src/extprime/pipelines/pipeline_supervised.py \
  --exp_type exp1 --max_k 6 --data_root data --output_root outputs
```

### 2) Accuracy vs. Genome Proxy Length (Single Run)

Compares accuracy across proxy lengths (k set by `--max_k`).

```bash
python3 src/extprime/pipelines/pipeline_supervised.py \
  --exp_type exp2 --max_k 6 --data_root data --output_root outputs
```

### 3) Effect of Number of Subfragments (n)

Varies n in the composite genome proxy to measure its impact.

```bash
python3 src/extprime/pipelines/pipeline_supervised.py \
  --exp_type exp3 --max_k 6 --data_root data --output_root outputs
```

**Optional:** Add `--whole_genome` to use entire genomes instead of proxies.

## Outputs

Results are written under:

```
outputs/{exp_type}/{env}/fragments_L{length}[_n{n}]/...
```

Each folder contains the generated FASTA for that environment and the model outputs produced by the pipeline.

## Command-Line Flags

- `--exp_type {exp1,exp2,exp3,tuning}` – Choose the experiment
- `--max_k INT` – Maximum k-mer length considered by the models
- `--data_root PATH` – Input root (default: `data`)
- `--output_root PATH` – Results root (default: `outputs`)
- `--whole_genome` – Use entire genomes instead of proxies (optional)

## Unupervised Learning Experiments

### 4) Non-parametric clustering and candidate identification

```bash
python3 src/extprime/pipelines/pipeline_unsupervised.py \
  --exp_type non-parametric --k_mer 6 --data_root path_to_the_subfragments --output_root outputs \
  --fragement_length 100000 --n_clusters 4 --env Temperature
```
## Command-Line Flags
- `--exp_type {parametric, non-parametric} – Choose the experiment`
- `--k_mer INT – K-mer length`
- `--fragement_length INT – Fragment length`
- `--n_clusters INT – Number of clusters (default: 4) - not needed for non-parametric`
- `--outputs_root PATH – Output results directory`
- `--env {pH,Temperature} – Environment type`
- `--fragment_path PATH – Path to a fragment FASTA file`
- `--data_root PATH – Input data directory`

## FCGR distance calculation and filtering

```bash
python3 src/extprime/analysis/distance_calculator.py --data_root path_to_the_subfragments
```



