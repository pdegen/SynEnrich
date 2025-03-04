# scripts/combine_libs.py

import os
import yaml
import pandas as pd
from typing import List, Dict
from explore_results import get_sig_dict, create_intersection_depth_df
from utils import pickler, read_enrichr, read_gmt


def create_summary_dict(
    savepath: str,
    lib_names: Dict[str, str],
    tools: List[str],
    metrics: List[str],
    project_name: str,
    qval: float = 0.05,
    save: bool = False,
    go_sem_sim: bool = False,
) -> Dict:
    libs = lib_names.keys()
    summary_dict: Dict = {lib: dict() for lib in libs}
    for lib in libs:
        summary_df = pd.read_csv(f"{savepath}/syn.combined.{lib}.{project_name}.csv", index_col=0, header=[0, 1, 2])
        summary_df.sort_values(by=("Combined", "nan", "Combined FDR"))
        summary_dict[lib]["summary_df"] = summary_df
        sig_dict = get_sig_dict(summary_df, tools, metrics, qval=qval, verbose=True)
        depth_df = create_intersection_depth_df(sig_dict)
        depth_df = format_depth_df(depth_df, project_name, savepath, lib, lib_names, summary_df, go_sem_sim)
        summary_dict[lib]["depth_df"] = depth_df

    # Store results in dictionary for meta-analysis
    if save:
        print(f"Saving summary dict with qval threshold = {qval}")
        outfile = os.path.join(savepath, f"syn.summary_dict.{project_name}.txt")
        pickler(summary_dict, outfile)

    return summary_dict


def format_depth_df(
    depth_df: pd.DataFrame,
    project_name: str,
    savepath: str,
    lib: str,
    lib_names: Dict[str, str],
    summary_df: pd.DataFrame,
    go_sem_sim: bool = False,
) -> pd.DataFrame:
    outfile = os.path.join(savepath, f"syn.depth.{lib}.{project_name}.csv")
    d = depth_df
    if len(d) < 1:
        print(f"No terms found for {lib}")
        pd.DataFrame([f"No terms found for {lib}"]).to_csv(outfile)
        return d

    isGO = d.index[0].startswith("GO:")

    # Append GO sub-ontology
    if isGO:
        d["ID"] = d.index.str.split("_").str[0]
        ont = pd.DataFrame(summary_df[("nan", "nan", "ONTOLOGY")].values, index=summary_df.index, columns=["ONTOLOGY"])
        d = d.merge(ont["ONTOLOGY"], left_on="ID", right_index=True, how="left")

    if "Direction" not in d:
        d["Direction"] = d.index.str.split("_").str[1].str.strip()

    d.index = d.index.str.split("_").str[0]
    d.rename({"Factors": "Configurations"}, axis=1, inplace=True)
    d["Combined FDR"] = summary_df.loc[d.index, ("Combined", "nan", "Combined FDR")]

    d.sort_values(by=["Depth", "Combined FDR"], ascending=[False, True], inplace=True)

    # Column order
    cols = ["Description", "Depth", "Direction", "Combined FDR"]

    if "ONTOLOGY" in d.columns:
        cols.append("ONTOLOGY")

    # Check which terms are in Enrichr library
    enrichr = "resources/Ontologies/GO_Enrichr_2023.gmt"  # TO DO: pass as arg
    if os.path.isfile(enrichr):
        print("Enrichr gmt file found, adding info to depth df")
        enrichr_df = read_enrichr(enrichr)
        d["Enrichr"] = d.index.isin(enrichr_df.index)
        cols.append("Enrichr")

    # Needs different conda env!
    # Check if we should do GO Semantic Similarity Analysis
    # if go_sem_sim and isGO:
    #     from clustering import append_GO_clusters_to_depth_df
    #     d = append_GO_clusters_to_depth_df(sim_matrix_cache_folder, d, figpath, lib, project_name, max_thresh=go_sem_sim_max_distance)
    #     cols.append("GO_Cluster")
    #     cols.append("Top_GO_Cluster")

    cols.append("Configurations")

    # Check if genes from gmt file can be appended
    gmt_file = lib_names[lib]
    gmt_found = False
    if gmt_file.endswith(".gmt") and os.path.isfile(gmt_file):
        gmt_found = True
    elif os.path.isfile(os.path.join("resources/Ontologies", gmt_file)):
        gmt_file = os.path.join("resources/Ontologies", gmt_file)
        gmt_found = True
    if gmt_found:
        print("Attempting to append genes from .gmt to depth df")
        gmt = read_gmt(gmt_file)
        if "Genes" in gmt:
            gmt.set_index("ID", inplace=True)
            d = d.merge(gmt["Genes"], left_index=True, right_index=True, how="left")
            d["Genes"] = d["Genes"].apply(lambda s: "; ".join(g for g in s))
            cols.append("Genes")
        else:
            print("No 'Genes' column found in gmt file...")

    d = d[cols]
    d.to_csv(outfile)
    return d


if __name__ == "__main__":
    config_file_path = os.path.join("config", "config.yaml")
    with open(config_file_path, "r") as stream:
        config = yaml.safe_load(stream)

    libs = config.get("libraries", [])
    lib_names = {lib.split(".gmt")[0]: lib for lib in libs}
    tools = config.get("tools", [])
    metrics = config.get("metrics", [])
    project_name = config.get("project_name")
    qval = config.get("qval")
    save = config.get("save_summary_dict")
    savepath = os.path.join("results", project_name, "combined")

    # GoSemSim
    go_sem_sim = config.get("go_sem_sim")
    if go_sem_sim:
        figpath = os.path.join("results", project_name, "figures")
        organismKEGG = config.get("organismKEGG")
        go_sem_sim_max_distance = config.get("go_sem_sim_max_distance")
        sim_matrix_cache_folder = os.path.join("results", project_name, ".cache")

        match organismKEGG:
            case "mmu":
                org = "Org.Mm.eg.db"
            case "hsa":
                org = "Org.Hs.eg.db"
            case _:
                raise Exception(f"Organism not supported: {organismKEGG}")

    create_summary_dict(savepath, lib_names, tools, metrics, project_name, qval, save, go_sem_sim=go_sem_sim)
