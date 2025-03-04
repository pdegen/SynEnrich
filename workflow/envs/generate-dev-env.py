# Annoying workaround since conda doesn't support including other environment.yaml files, unlike --requirement requirements.txt

import yaml


# Define the dev packages to append
DEV_PACKAGES = [
    #   "setuptools>=68.2.2",
    #   "build>=1.0.3",
    #   "twine>=4.0.2",
    # Linting/Tooling
    # "black>=23.11.0",
    "isort>=5.12.0",
    # "mypy>=1.6.1",
    "pre-commit>=3.5.0",
    "ruff>=0.1.5",
    "pyright>=1.1.335",
    # Testing
    "pytest>=7.4.3",
    "codecov>=2.1.13",
    "pytest-cov>=4.1.0",
    "pytest-benchmark>=4.0.0",
    # Documentation
    "mkdocs>=1.5.3",
    "mkdocstrings>=0.23.0",
    "mkdocstrings-python>=1.7.3",
    "mkdocs-material>=9.4.8",
    "Pygments>=2.16.1",

    # Project-specific
    "snakemake=8.18.1",
    "mygene"
]

# Load the base environment.yaml
with open("environment.yaml", "r") as f:
    base_env = yaml.safe_load(f)

# Create a copy for the dev environment
dev_env = base_env.copy()
dev_env["name"] = "SynDev"

# Append dev dependencies
if "dependencies" not in dev_env:
    dev_env["dependencies"] = []

dev_env["dependencies"].extend(DEV_PACKAGES)

# Save the new environment-dev.yaml
with open("environment.dev.yaml", "w") as f:
    yaml.dump(dev_env, f, default_flow_style=False)

print("✅ Generated environment.dev.yaml successfully!")
