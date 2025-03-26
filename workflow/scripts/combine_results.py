# scripts/combine_results.py

import os
import sys
import argparse
import glob
from typing import Optional, List, Union
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from scripts.utils import load_config


def stouffer_combined_p_value(
    p_values: Union[List[float], np.ndarray], weights: Optional[Union[List[float], np.ndarray]] = None
) -> np.ndarray:
    if weights is None:
        weights = [1] * len(p_values)
    else:
        weights = np.array(weights)

    z_enrichmentScores = norm.ppf(1 - np.array(p_values))
    weighted_z = np.sum(weights * z_enrichmentScores) / np.sqrt(np.sum(np.array(weights) ** 2))
    combined_p_value = norm.cdf(-weighted_z)
    return combined_p_value


def mean_pval(pvals):
    pv = pvals.copy()  # pvals is read-only
    pv[pv == 0] = 1e-30  # Impute zeros with a small number; TO DO: change to 1/n_perms or something
    pv = pv.astype(float)  # else np.log10 might complain
    return 10 ** np.nanmean(np.log10(pv))


## TO DO: refactor
def combine_results(dfs, isGO: bool = False) -> pd.DataFrame:
    combined_df = pd.concat(dfs, axis=1, keys=[d.index.name for d in dfs])
    combined_df.columns = pd.MultiIndex.from_product([combined_df.columns.levels[0], combined_df.columns.levels[1]])

    # Calculate combined pvals only for common terms
    combined_df_common = combined_df  # .dropna(axis=0)
    summary_df = pd.DataFrame(index=combined_df_common.index)
    summary_df["enrichmentScore Mean"] = combined_df_common.xs("enrichmentScore", axis=1, level=1).apply(
        np.nanmean, axis=1
    )
    summary_df["enrichmentScore SD"] = combined_df_common.xs("enrichmentScore", axis=1, level=1).apply(
        np.nanstd, axis=1
    )
    summary_df["Combined pvalue"] = combined_df_common.xs("pvalue", axis=1, level=1).apply(mean_pval, axis=1)
    summary_df["Combined FDR"] = fdrcorrection(
        summary_df["Combined pvalue"],
    )[1]
    summary_df.columns = pd.MultiIndex.from_product([["Combined"], summary_df.columns])
    summary_df = pd.concat([combined_df_common, summary_df], axis=1)

    # Add remaining terms that are not common to all tools
    # combined_df_uncommon = combined_df.loc[combined_df.index.difference(combined_df_common.index)]
    # summary_df_uncommon = pd.DataFrame(index=combined_df_uncommon.index)
    # summary_df_uncommon["enrichmentScore Mean"] = combined_df_uncommon.xs("enrichmentScore", axis=1, level=1).mean(axis=1)
    # summary_df_uncommon["enrichmentScore SD"] = combined_df_uncommon.xs("enrichmentScore", axis=1, level=1).std(axis=1)
    # summary_df_uncommon["Combined pvalue"] = np.nan
    # summary_df_uncommon["Combined pvalue"] = np.nan
    # summary_df_uncommon.columns = pd.MultiIndex.from_product([["Combined"], summary_df_uncommon.columns])
    # summary_df_uncommon = pd.concat([combined_df_uncommon,summary_df_uncommon], axis=1)

    # summary_df = pd.concat([summary_df, summary_df_uncommon], axis=0)

    def combine_col(row):
        vals = row.dropna().unique()
        return vals[0] if vals.size > 0 else np.nan

    summary_df["Description"] = summary_df.loc[:, pd.IndexSlice[:, "Description"]].apply(combine_col, axis=1)
    summary_df.drop("Description", level=1, axis=1, inplace=True)

    if isGO:
        summary_df["ONTOLOGY"] = summary_df.loc[:, pd.IndexSlice[:, "ONTOLOGY"]].apply(combine_col, axis=1)
        summary_df.drop("ONTOLOGY", level=1, axis=1, inplace=True)

    # Better multiindex
    lvl0 = summary_df.columns.get_level_values(level=0)
    lvl0 = [val.split(".")[0] for val in lvl0]
    lvl00 = summary_df.columns.get_level_values(level=0)
    lvl1 = [val.split(".")[1] if len(val.split(".")) > 1 else "nan" for val in lvl00]
    lvl2 = list(summary_df.columns.get_level_values(level=1))

    lvl0[-1] = "nan"

    if isGO:
        lvl0[-2] = "nan"
        lvl2[-2] = "Description"
        lvl2[-1] = "ONTOLOGY"
    else:
        lvl2[-1] = "Description"

    tuples = list(zip(*[lvl0, lvl1, lvl2], strict=True))
    summary_df.columns = pd.MultiIndex.from_tuples(tuples, names=["Tool", "Metric", "Value"])

    return summary_df


def format_table(tab: pd.DataFrame, tool: str, metric: str, library: str) -> pd.DataFrame:
    tab.rename(
        {"NOM p-val": "pvalue", "ES": "enrichmentScore", "Term": "Description", "FDR q-val": "qvalue"},
        axis=1,
        inplace=True,
    )

    # currently unused
    # tab["Direction"] = tab["enrichmentScore"].apply(lambda x: "Up" if x > 0 else "Down")
    # tab["Signed_Term"] = tab.index + "." + tab["Direction"]

    if library == "KEGG":
        # TO DO: fix this in clusterprofiler script
        if tool == "clusterProfiler" and " - Mus musculus (house mouse)" in tab["Description"].iloc[0]:
            tab["Description"] = tab["Description"].str.replace(" - Mus musculus (house mouse)", "")
        tab.set_index("Description", drop=False, inplace=True)  # gseapy doesn't store KEGG IDs...

    elif tab.index.name != "ID":
        tab.set_index("ID", drop=False)

    tab.index.name = tool + "." + metric  # hacky
    return tab


def main(savepath: str, output_files: List[str], project_name: str) -> None:
    # Get latest config file
    config_file = os.path.join(savepath, "config.yaml")
    orig_config_file = os.path.join("config", "config.yaml")
    os.system(f"cp {orig_config_file} {config_file}")
    config = load_config(config_file)

    # metrics = config.get("metrics", [])
    libraries = config.get("libraries", [])
    tools = config.get("tools", [])

    input_files = glob.glob(f"{savepath}/syn.*[tc]sv")

    if len(input_files) < 2:
        print("Savepath:", savepath)
        print("Fewer than 2 input_files found, nothing to combine...")
        return
    else:
        print(f"Found {len(input_files)} input files:\n", *[o + "\n" for o in input_files])

    for library in libraries:
        if library.endswith(".gmt"):
            library = library.split(".gmt")[0]

        tab_dict = {}
        output_files_lib = next(o for o in output_files if library in o)  # TO DO: careful

        isGO = library.startswith("GO")  # TO DO: careful

        for tool in tools:
            for file in input_files:
                filename = os.path.basename(file)
                file_tool, file_metric, file_lib = filename.split("syn.")[1].split(".")[:3]

                if library != file_lib or tool != file_tool:
                    continue

                metric = file_metric

                sep = "\t" if os.path.splitext(file)[-1] == ".tsv" else ","
                tab = pd.read_csv(file, index_col=0, sep=sep)
                tab = format_table(tab, tool, metric, library)
                tab_dict[tab.index.name] = tab

                # string contains terms enriched in "both" directions;
                # also top refers to negative value rankings and vice versa
                if tool.lower() == "string":
                    tab["Direction"] = tab["direction"].apply(
                        lambda x: "Up" if x == "bottom" else "Down" if x == "top" else "Both"
                    )
                else:
                    tab["Direction"] = tab["enrichmentScore"].apply(lambda x: "Up" if x > 0 else "Down")

                if isGO and "Ontology" in tab:
                    tab.rename({"Ontology": "ONTOLOGY"}, axis=1, inplace=True)

        # Combine results
        cols = ["enrichmentScore", "pvalue", "qvalue", "Description", "Direction"]
        if isGO:
            cols += ["ONTOLOGY"]

        dfs = [d[cols] for d in tab_dict.values()]
        summary_df = combine_results(dfs, isGO)
        summary_df.to_csv(output_files_lib, index=True)


if __name__ == "__main__":
    pd.options.mode.copy_on_write = True

    parser = argparse.ArgumentParser(description="Combine results and save output.")

    # Positional argument for the save path
    parser.add_argument("savepath", type=str, help="Path to save the combined results.")

    # Positional argument for a variable number of final outputs
    parser.add_argument("final_outputs", nargs="+", help="List of final output files.")

    args = parser.parse_args()
    savepath = args.savepath
    output_files = args.final_outputs
    project_name = savepath.split("results/")[1]

    main(savepath, output_files, project_name)
