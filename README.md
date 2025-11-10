# MAG Functional Annotation Pipeline

A Snakemake workflow for bacterial and archaeal MAG functional annotation.

## Features
- Predicts genes using **Prokka**
- Annotates proteins with **eggNOG-mapper**
- Summarizes **KEGG KOs** and **COG categories**
- Generates combined KO and COG abundance tables

## Quick Start
```bash
conda install -c bioconda snakemake
git clone https://github.com/<your_username>/mag_functional_pipeline.git
cd mag_functional_pipeline
snakemake --use-conda --cores 8
