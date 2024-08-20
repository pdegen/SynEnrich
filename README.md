# SynEnrich

This package provides a convenient Snakemake workflow to synthesize functional scoring enrichment results from several sources (GSEApy, ClusterProfiler, STRING) and ranking metrics (logFC, signed p-value).

## Introduction

Practitioners of functional scoring analysis (commonly known as gene set enrichment analysis or GSEA) face two main challenges: 1) a proliferation of enrichment tools, which makes it difficult to determine the best approach, and 2) a lack of consensus on the optimal ranking metrics. Given these issues, a sensible strategy is to employ multiple methods and metrics and integrate the results to achieve more robust results.

## Usage

The notebook [main/workflow/notebooks/main.ipynb](https://github.com/pdegen/SynEnrich/blob/main/workflow/notebooks/main.ipynb) serves as the entry point. Users define variables (input files, parameters) in the notebook, which generates a config.yaml file for Snakemake. Then, the workflow must be run from the project root directory: 

`snakemake --use-conda --cores 1` (adjust number of cores as needed)

The results are stored as a csv file and can be inspected in the same notebook, which provides dedicated plotting methods.
