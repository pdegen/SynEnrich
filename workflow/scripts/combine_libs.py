# scripts/combine_results.py

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
        summary_dict[lib]["depth_df"] = create_intersection_depth_df(sig_dict)    

    # Store results in dictionary for meta-analysis
    if save:
        print(f"Saving summary dict with qval threshold = {qval}")
        outfile = os.path.join(savepath, f"syn.summary_dict.{project_name}.txt")
        pickler(summary_dict, outfile)

    return summary_dict

if __name__ == "__main__":

    config_file_path = os.path.join("config", "config.yaml")
    with open(config_file_path, 'r') as stream:
        config = yaml.safe_load(stream)

    libs = config.get('libraries', [])
    tools = config.get('tools', [])
    metrics = config.get('metrics', [])
    project_name = config.get('project_name')
    qval = config.get('qval')
    save = config.get('save_summary_dict')
    savepath = os.path.join("results",project_name,"combined")

    create_summary_dict(savepath, libs, tools, metrics, project_name, qval, save)



