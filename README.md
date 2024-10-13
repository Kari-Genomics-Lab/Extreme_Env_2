Note: The repo is not finalized yet. 
# Extremophiles Genomic Signature Analysis

This repository contains the code and data for the paper titled "Maximally divergent microbes share similar genomic signatures linked to extreme environments" by Monireh Safari, Joseph Butler, Gurjit S Randhawa, Kathleen A Hill, and Lila Kari.

## Overview

This study uncovers new microbial organism pairs whose genomic signatures are similar, despite their maximal taxonomic differences. The research introduces a pipeline for optimizing DNA representative fragment selection and k-mer length determination in genomic signature analysis of extremophiles. By employing an approach to select representative DNA fragments, the aim is to achieve a randomized and unbiased genomic representation essential for accurate taxonomy and environmental signature deciphering.

## Key Contributions

1. **Optimized DNA Fragment Selection:** 
   - Developed a method to select representative DNA fragments that ensures randomness and unbiased representation.
2. **k-mer Length and Fragment Size Optimization:** 
   - Determined the optimal k-mer length and DNA fragment size for classification accuracy and computational efficiency.
3. **Genomic Signature Pervasiveness:** 
   - Verified that genomic signatures are pervasive throughout the genome, ensuring high accuracy regardless of the selected genomic region.
4. **Microbial Pairs Identification:** 
   - Identified microbial species pairs with similar genomic signatures related to extreme environments despite their taxonomic differences.

## Repository Structure

```
├── data/
│   ├── raw/                  # Raw data files
│   └── processed/            # Processed data files
├── results/                  # Results from the analysis
├── scripts/
│   ├── analysis/
│   │   ├── analyse_results_exp2.py
│   │   ├── candidates_dist_analysis.py
│   │   └── distance_analysis.py
│   ├── download/
│   │   └── download_assemblies.py
│   ├── models/
│   │   ├── supervised_models.py
│   │   ├── supervised_models_challenging.py
│   │   ├── unsupervised_non_parametric.py
│   │   └── unsupervised_parametric.py
│   ├── pipelines/
│   │   ├── pipeline_supervised.py
│   │   └── pipeline_unsupervised.py
└── README.md                 # This README file
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/msafari18/Ext_env_2.git
   cd Ext_env_2
   ```

2. Create and activate a virtual environment:
   ```bash
   python3 -m venv env
   source env/bin/activate
   ```

3. Install the required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Data Preparation

1. Download the required genome assemblies using the script in the `download` directory:
   ```bash
   python scripts/download/download_assemblies.py
   ```

2. Build the signature dataset:
   ```bash
   python scripts/build_signature_dataset.py
   ```

### Analysis

1. Run the supervised pipeline:
   ```bash
   python scripts/pipelines/pipeline_supervised.py
   ```

2. Run the unsupervised pipeline:
   ```bash
   python scripts/pipelines/pipeline_unsupervised.py
   ```

3. Perform distance analysis:
   ```bash
   python scripts/analysis/distance_analysis.py
   ```

4. Analyze the results:
   ```bash
   python scripts/analysis/analyse_results_exp2.py
   ```

## Results

The results of the analysis, including the optimized k-mer lengths, DNA fragment sizes, and identified microbial pairs, are stored in the `results` directory.

## Figures

### Figure 3
![Figure 3: Multi-Layered Filter for identifying bacteria/archaea pairs with similar genomic signatures.](fig3.png)

### Figure 4
![Figure 4: FCGR images of 2 pairs (4 species) from the final set of 11 pairs, with a resolution of k = 8.](fig4.png)

## Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada Discovery Grants and Compute Canada Research Platforms & Portals Grant.

## Contact

For any questions or issues, please contact Monireh Safari at monireh.safari@uwaterloo.ca.
