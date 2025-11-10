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
git clone https://github.com/bhimbbiswa/MAG_functional_annotation.git
cd mag_functional_pipeline
snakemake --use-conda --cores 8
