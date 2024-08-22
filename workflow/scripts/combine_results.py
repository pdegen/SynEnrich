# scripts/combine_results.py

import os
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
from matplotlib_venn import venn2, venn3
import glob
from typing import Optional, List, Union, Dict, Any

def stouffer_combined_p_value(p_values: Union[List[float], np.ndarray], 
                              weights: Optional[Union[List[float], np.ndarray]] = None) -> float:
    if weights is None:
        weights = [1] * len(p_values)
    else:
        weights = np.array(weights)
        
    z_enrichmentScores = norm.ppf(1 - np.array(p_values))
    weighted_z = np.sum(weights * z_enrichmentScores) / np.sqrt(np.sum(np.array(weights) ** 2))
    combined_p_value = norm.cdf(-weighted_z)
    return combined_p_value

def combine_results(dfs: List[pd.DataFrame]) -> pd.DataFrame:

    combined_df = pd.concat(dfs, axis=1, keys=[d.index.name for d in dfs])
    combined_df.columns = pd.MultiIndex.from_product([combined_df.columns.levels[0], combined_df.columns.levels[1]])

    summary_df = pd.DataFrame(index=combined_df.index)
    summary_df["enrichmentScore Mean"] = combined_df.xs("enrichmentScore", axis=1, level=1).mean(axis=1)
    summary_df["enrichmentScore SD"] = combined_df.xs("enrichmentScore", axis=1, level=1).std(axis=1)
    summary_df["Stouffer pvalue"] = combined_df.xs("pvalue", axis=1, level=1).apply(stouffer_combined_p_value, axis=1)
    summary_df["Stouffer FDR"] = fdrcorrection(summary_df["Stouffer pvalue"], )[1]

    summary_df.columns = pd.MultiIndex.from_product([["Combined"], summary_df.columns])

    return pd.concat([combined_df,summary_df], axis=1)

def load_config(file_path: str) -> Dict[str, Any]:
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main(savepath, output_files, project_name):

    config_file = f"{savepath}/config.yaml"
    try:
        config = load_config(config_file)
    except FileNotFoundError:
        print("No local config file found, trying to copy from root/config")
        config = load_config("config/config.yaml")
        config_project_name = config.get('project_name')
        
        if config_project_name == project_name:
            ### Copy config file
            os.system(f"cp config/config.yaml {config_file}")
        else:
            raise Exception("Root config doesn't match project_name; re-generate config.yaml")

    metrics = config.get('metrics', [])
    libraries = config.get('libraries', [])
    tools = config.get('tools', [])

    input_files = glob.glob(f"{savepath}/syn.*csv")

    if len(input_files) < 2:
        print("Savepath:", savepath)
        print("Fewer than 2 input_files found, nothing to combine...")
        return
    else:
        print(f"Found {len(input_files)} input files:\n",*[o+"\n" for o in input_files])

    for library in libraries:

        tab_dict = dict()

        output_files_lib = [o for o in  output_files if library in o][0] # TO DO: careful

        for tool in tools:
            
            for file in input_files:

                if library not in file or tool not in file: 
                    continue

                metric = [m for m in metrics if m in file][0]  # TO DO: careful
                tab = pd.read_csv(file, index_col=0)
                tab["Direction"] = tab["enrichmentScore"].apply(lambda x: "Up" if x > 0 else "Down")
                tab["Signed_Term"] = tab.index + "_" + tab["Direction"]
                tab.index.name = tool + "_" + metric # hacky
                tab_dict[tab.index.name] = tab
        
        ### Combine results
        dfs = [d[["enrichmentScore","pvalue"]] for d in tab_dict.values()]
        summary_df = combine_results(dfs)

        summary_df["Description"] = tab_dict[tab.index.name].loc[summary_df.index,"Description"]
        summary_df.to_csv(output_files_lib, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine results and save output.")
    
    # Positional argument for the save path
    parser.add_argument("savepath", type=str, help="Path to save the combined results.")
    
    # Positional argument for a variable number of final outputs
    parser.add_argument("final_outputs", nargs='+', help="List of final output files.")
    
    args = parser.parse_args()
    savepath = args.savepath
    output_files = args.final_outputs
    project_name = savepath.split("results/")[1]

    main(savepath, output_files, project_name)



