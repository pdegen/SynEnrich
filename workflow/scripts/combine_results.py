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

def combine_results2(dfs: List[pd.DataFrame]) -> pd.DataFrame:

    combined_df = pd.concat(dfs, axis=1, keys=[d.index.name for d in dfs])
    combined_df.columns = pd.MultiIndex.from_product([combined_df.columns.levels[0], combined_df.columns.levels[1]])

    print("Before:", len(combined_df))
    combined_df.dropna(axis=0, inplace=True)
    print("After:", len(combined_df))

    summary_df = pd.DataFrame(index=combined_df.index)
    summary_df["enrichmentScore Mean"] = combined_df.xs("enrichmentScore", axis=1, level=1).mean(axis=1)
    summary_df["enrichmentScore SD"] = combined_df.xs("enrichmentScore", axis=1, level=1).std(axis=1)
    summary_df["Stouffer pvalue"] = combined_df.xs("pvalue", axis=1, level=1).apply(stouffer_combined_p_value, axis=1)
    summary_df["Stouffer FDR"] = fdrcorrection(summary_df["Stouffer pvalue"], )[1]

    summary_df.columns = pd.MultiIndex.from_product([["Combined"], summary_df.columns])

    return pd.concat([combined_df,summary_df], axis=1)

## TO DO: refactor
def combine_results(dfs) -> pd.DataFrame:

    combined_df = pd.concat(dfs, axis=1, keys=[d.index.name for d in dfs])
    combined_df.columns = pd.MultiIndex.from_product([combined_df.columns.levels[0], combined_df.columns.levels[1]])

    # Calculate combined pvals only for common terms
    combined_df_common = combined_df.dropna(axis=0)
    summary_df = pd.DataFrame(index=combined_df_common.index)
    summary_df["enrichmentScore Mean"] = combined_df_common.xs("enrichmentScore", axis=1, level=1).mean(axis=1)
    summary_df["enrichmentScore SD"] = combined_df_common.xs("enrichmentScore", axis=1, level=1).std(axis=1)
    summary_df["Stouffer pvalue"] = combined_df_common.xs("pvalue", axis=1, level=1).apply(stouffer_combined_p_value, axis=1)
    summary_df["Stouffer FDR"] = fdrcorrection(summary_df["Stouffer pvalue"], )[1]
    summary_df.columns = pd.MultiIndex.from_product([["Combined"], summary_df.columns])
    summary_df = pd.concat([combined_df_common,summary_df], axis=1)

    # Add remaining terms that are not common to all tools
    combined_df_uncommon = combined_df.loc[combined_df.index.difference(combined_df_common.index)]
    summary_df_uncommon = pd.DataFrame(index=combined_df_uncommon.index)
    summary_df_uncommon["enrichmentScore Mean"] = combined_df_uncommon.xs("enrichmentScore", axis=1, level=1).mean(axis=1)
    summary_df_uncommon["enrichmentScore SD"] = combined_df_uncommon.xs("enrichmentScore", axis=1, level=1).std(axis=1)
    summary_df_uncommon["Stouffer pvalue"] = np.nan
    summary_df_uncommon["Stouffer pvalue"] = np.nan
    summary_df_uncommon.columns = pd.MultiIndex.from_product([["Combined"], summary_df_uncommon.columns])
    summary_df_uncommon = pd.concat([combined_df_uncommon,summary_df_uncommon], axis=1)
    
    summary_df = pd.concat([summary_df, summary_df_uncommon], axis=0)

    def combine_descriptions(row):
        descriptions = row.dropna().unique()
        return ', '.join(descriptions) if descriptions.size > 0 else np.nan

    summary_df['Description'] = summary_df.loc[:, pd.IndexSlice[:, 'Description']].apply(combine_descriptions, axis=1)
    summary_df.drop("Description", level=1, axis=1, inplace=True)
    
    # Better multiindex
    lvl0 = summary_df.columns.get_level_values(level = 0)
    lvl0 = [val.split(".")[0] for val in lvl0]
    lvl00 = summary_df.columns.get_level_values(level = 0)
    lvl1 = [val.split(".")[1] if len(val.split("."))>1 else "" for val in lvl00]
    lvl2 = list(summary_df.columns.get_level_values(level = 1))
    lvl0[-1] = ""
    lvl2[-1] = "Description"

    tuples = list(zip(*[lvl0,lvl1,lvl2]))
    summary_df.columns = pd.MultiIndex.from_tuples(tuples, names=["Tool","Metric","Value"])

    return summary_df

def load_config(file_path: str) -> Dict[str, Any]:
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main(savepath: str, output_files: List[str], project_name: str) -> None:

    # Get latest config file
    config_file = f"{savepath}/config.yaml"
    os.system(f"cp config/config.yaml {config_file}")
    config = load_config(config_file)


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
                tab["Signed_Term"] = tab.index + "." + tab["Direction"]

                if library == "GO":
                    if tab.index.name != "ID":
                        tab.index = tab["ID"]
                elif library == "KEGG":
                    tab.index = tab["Description"] # gseapy doesn't store KEGG IDs...
                else:
                    raise Exception("Library not supported:", library)
                    
                tab.index.name = tool + "." + metric # hacky
                tab_dict[tab.index.name] = tab
        
        ### Combine results
        dfs = [d[["enrichmentScore","pvalue","Description"]] for d in tab_dict.values()]
        summary_df = combine_results(dfs)
        summary_df.to_csv(output_files_lib, index=True)

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



