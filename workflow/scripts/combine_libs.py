# scripts/combine_libs.py

import os
import yaml
import pandas as pd
from typing import List, Dict
from explore_results import get_sig_dict, create_intersection_depth_df
from utils import pickler

def create_summary_dict(savepath: str,
         libs: List[str],
         tools: List[str],
         metrics: List[str],
         project_name: str,
         qval: float = 0.05,
         save = False,
         ) -> Dict:

    summary_dict: Dict = {lib: dict() for lib in libs}
    for lib in libs:
        summary_df = pd.read_csv(f"{savepath}/syn.combined.{lib}.{project_name}.csv", index_col=0, header=[0,1,2])
        summary_df.sort_values(by=("Combined","nan","Combined FDR"))
        summary_dict[lib]["summary_df"] = summary_df
        sig_dict = get_sig_dict(summary_df, tools, metrics, qval=qval, verbose=True)
        depth_df = create_intersection_depth_df(sig_dict)
        depth_df = format_depth_df(depth_df, project_name, savepath, lib, summary_df)
        summary_dict[lib]["depth_df"] = depth_df

    # Store results in dictionary for meta-analysis
    if save:
        print(f"Saving summary dict with qval threshold = {qval}")
        outfile = os.path.join(savepath, f"syn.summary_dict.{project_name}.txt")
        pickler(summary_dict, outfile)

    return summary_dict

def format_depth_df(depth_df: pd.DataFrame,
                    project_name: str,
                    savepath: str,
                    lib: str,
                    summary_df: pd.DataFrame
                    ) -> pd.DataFrame:
    
    outfile = os.path.join(savepath, f"syn.depth.{lib}.{project_name}.csv")
    d = depth_df
    if len(d) < 1:
        print(f"No terms found for {lib}")
        pd.DataFrame(f"No terms found for {lib}").to_csv(outfile)
        return d
    
    # Append GO sub-ontology
    if d.index[0].startswith("GO:"):
        d["ID"] = d.index.str.split("_").str[0]
        ont = pd.DataFrame(summary_df[("nan","nan","ONTOLOGY")].values, index=summary_df.index, columns=["ONTOLOGY"])
        d = d.merge(ont["ONTOLOGY"], left_on="ID", right_index=True, how="left")

    if "Direction" not in d:
        d["Direction"] = d.index.str.split("_").str[1].str.strip()

    d.index = d.index.str.split("_").str[0]
    d.rename({"Factors":"Configurations"}, axis=1, inplace=True)
    d["Combined FDR"] = summary_df.loc[d.index,("Combined","nan","Combined FDR")]

    d.sort_values(by=["Depth","Combined FDR"], ascending=[False,True], inplace=True)

    cols = ["Description","Depth","Direction","Combined FDR"]
    if "ONTOLOGY" in d.columns:
        cols.append("ONTOLOGY")
    cols.append("Configurations")

    d = d[cols]
    d.to_csv(outfile)
    return d

if __name__ == "__main__":

    config_file_path = os.path.join("config", "config.yaml")
    with open(config_file_path, 'r') as stream:
        config = yaml.safe_load(stream)

    libs = config.get('libraries', [])
    libs = [lib.split(".gmt")[0] for lib in libs]
    tools = config.get('tools', [])
    metrics = config.get('metrics', [])
    project_name = config.get('project_name')
    qval = config.get('qval')
    save = config.get('save_summary_dict')
    savepath = os.path.join("results",project_name,"combined")

    create_summary_dict(savepath, libs, tools, metrics, project_name, qval, save)




