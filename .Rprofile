#### note: when uncommenting this, jobs don't use correct .snakemake conda envs

# # .Rprofile
# # Dynamically set .libPaths based on the active Conda environment
# print("Loading .Rprofile...")

# # Retrieve the CONDA_PREFIX environment variable, which points to the active Conda environment
# conda_prefix <- Sys.getenv("CONDA_PREFIX")

# if (nzchar(conda_prefix)) {
#   # Construct the path to the R library within the Conda environment
#   conda_lib_path <- file.path(conda_prefix, "lib", "R", "library")

#   # Set .libPaths to use only the Conda environment's library path
#   .libPaths(conda_lib_path)
# } else {
#   warning("CONDA_PREFIX is not set. Using default .libPaths.")
# }

# # Optionally print the new .libPaths to confirm
# print(.libPaths())