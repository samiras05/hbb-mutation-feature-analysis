# HBB Thalassemia Mutation Portfolio

This repository contains a clean, thesis-inspired bioinformatics workflow for analyzing HBB missense mutations in beta-thalassemia.

The project focuses on:
- amino acid biochemical features
- Grantham and BLOSUM62 scoring
- physicochemical mutation shifts
- DynaMut stability outputs
- FoldX energy terms
- integration of structural and genomic predictors
- correlation analysis, PCA, and clustering

## Files

- `scripts/01_feature_analysis.py`
- `scripts/02_physicochemical_analysis.py`
- `scripts/03_dynamut_analysis.py`
- `scripts/04_foldx_analysis.py`
- `scripts/05_integrated_analysis.py`

## Input files

Place your input files in the project root or in a local `data/` folder:

- `HBB_155_merged_validated_ref_alt.xlsx`
- `AA.xlsx`
- `DynaMut.xlsx`
- `FoldX.xlsx`
- `MergedData.csv`

## Requirements

Install the Python packages listed in `requirements.txt`.

## How to run

```bash
python scripts/01_feature_analysis.py
python scripts/02_physicochemical_analysis.py
python scripts/03_dynamut_analysis.py
python scripts/04_foldx_analysis.py
python scripts/05_integrated_analysis.py
