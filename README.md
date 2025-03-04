# SynEnrich (pre-release version)

This package provides a convenient Snakemake workflow to synthesize functional scoring enrichment results from several sources (GSEApy, ClusterProfiler, STRING) and ranking metrics (e.g. logFC, signed p-value, signal-to-noise ratio, or any user-provided metric). A preprint of our study demonstrating this tool in more detailed will be published in the first half of 2025.

## Introduction

Practitioners of functional scoring analysis (commonly known as gene set enrichment analysis or GSEA) face two main challenges: 1) a proliferation of enrichment tools, which makes it difficult to determine the best method, and 2) a lack of consensus on the optimal ranking metric. Given these issues, a sensible strategy is to employ multiple tools and metrics to assess the robustness of results. We designate a unique combination of software tool and ranking metric as an analysis _configuration_. Reproducible workflow manager such as Snakemake provide a convenient and streamlined way to run multiple configurations in parallel, as well as to add more configurations and update the results dynamically.

## Usage

The notebook [main/workflow/notebooks/main.ipynb](https://github.com/pdegen/SynEnrich/blob/main/workflow/notebooks/main.ipynb) serves as the entry point. Users define variables (input files, parameters) in the notebook, which generates a config.yaml file for Snakemake. Then, the workflow should be run from the project root directory:

`snakemake --use-conda --cores 1` (adjust number of cores as needed)

The workflow generates a CSV file with all significant enrichment terms ranked by their _robustness_, i.e., the number of analysis configurations in which a given term is found significant. Additionally, various figures are generated to visualize the results (bar plots, Venn diagrams, UpSet plots).

## To do

- Perform STRING values/ranks enrichment using recent API update
