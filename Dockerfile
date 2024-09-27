FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2d9d6cbf72b9938a1bdd662db4090128b0f944668daf529aad5efaa706338f48"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/environment.clusterprofiler.yaml
#   prefix: /conda-envs/9e5d21de02f9a63b44a8ae2ec3a88b9f
#   # environment.yaml
#   
#   name: SynEnrich_ClusterProfiler
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - mamba
#     - r-base=4.3.3
#     - r-essentials
#     - bioconductor-clusterprofiler=4.10.0
#     - bioconductor-org.hs.eg.db=3.18.0
#     - bioconductor-org.mm.eg.db=3.18.0
#     - bioconductor-annotationdbi=1.64.1
RUN mkdir -p /conda-envs/9e5d21de02f9a63b44a8ae2ec3a88b9f
COPY workflow/envs/environment.clusterprofiler.yaml /conda-envs/9e5d21de02f9a63b44a8ae2ec3a88b9f/environment.yaml

# Conda environment:
#   source: workflow/envs/environment.yaml
#   prefix: /conda-envs/145ca9dde66d47e17a55541c7109c5f8
#   # environment.yaml
#   
#   name: SynEnrich
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - python=3.12.3
#     - mamba
#     - pip
#     - pyyaml
#     - pandas
#     - matplotlib
#     - seaborn
#     - statsmodels
#     - scipy
#     - matplotlib-venn
#     - upsetplot
#     - snakemake
#     - gseapy
#     - rpy2
#     - r-base=4.3.3
#     - r-essentials
#     - bioconductor-clusterprofiler=4.10.0
#     - bioconductor-org.hs.eg.db=3.18.0
#     - bioconductor-annotationdbi=1.64.1
RUN mkdir -p /conda-envs/145ca9dde66d47e17a55541c7109c5f8
COPY workflow/envs/environment.yaml /conda-envs/145ca9dde66d47e17a55541c7109c5f8/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/9e5d21de02f9a63b44a8ae2ec3a88b9f --file /conda-envs/9e5d21de02f9a63b44a8ae2ec3a88b9f/environment.yaml && \
    mamba env create --prefix /conda-envs/145ca9dde66d47e17a55541c7109c5f8 --file /conda-envs/145ca9dde66d47e17a55541c7109c5f8/environment.yaml && \
    mamba clean --all -y
