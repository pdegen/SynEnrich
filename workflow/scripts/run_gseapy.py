# scripts/run_gseapy.py

import os
import sys
import numpy as np
import pandas as pd
import gseapy
from typing import Union, List

def convert_gseapy_table(tab: pd.DataFrame, ont_id: str) -> None:
    '''Convert GSEApy columns to match clusterProfiler'''
    tab.rename({"NOM p-val": "pvalue", "ES": "enrichmentScore", 
                "FDR q-val": "qvalue"}, axis=1, inplace=True) # dubious combining scores?

    print("ont_id", ont_id)

    if ont_id == "KEGG":
        tab.rename({"Term": "Description"}, axis=1, inplace=True)
        tab["ID"] = tab["Description"] # gseapy doesn't save KEGG id, hence use description
    elif ont_id == "GO":
        print(tab.head())
        tab["ID"] = "GO:" + tab["Term"].str.split("\\(GO:").str[1].str[:-1]
    else:
        try:
            tab.rename({"Term":"ID"}, axis=1, inplace=True)
        except KeyError:
            pass
    
    if "ID" not in tab.columns:
        raise Exception("ID not found in results table, make sure gmt file has column called 'ID' or 'Term'")
    
    tab.set_index("ID", inplace=True, drop=False)

def run_gseapy_multi(   
    tab: pd.DataFrame, 
    metric: str, 
    ontology: str, # "GO", "KEGG", Enrichr library, or gmt file
    organismKEGG: str = "",
    outdir: str = "", 
    outfile: str = "",
    overwrite: bool = False,

    **kwargs) -> None:

    if metric not in tab.columns:
        if metric == "neg_signed_logpval":
            tab["neg_signed_logpval"] = -np.sign(tab["logFC"]) * np.log10(tab["PValue"])
        else:
            raise Exception(f"Metric not in input table: {metric}")
        
    input = tab[metric].sort_values(ascending=False)
    input.index.name = None # remove header

    print("Running:", ontology)

    if ontology == "GO":
        ont_id = "GO"
        ontology = ["GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"]
    elif ontology == "KEGG" and organismKEGG == "hsa":
        ont_id = "KEGG"
        ontology = "KEGG_2021_Human"
    elif ontology == "KEGG" and organismKEGG == "mmu":
        ont_id = "KEGG"
        ontology = "KEGG_2019_Mouse"
    elif ontology.endswith(".gmt"):
        ont_id = ontology.split("/")[-1].split(".gmt")[0]
        if not os.path.isfile(ontology):
            ontology = os.path.join("resources/Ontologies", ontology)
        if not os.path.isfile(ontology):
            raise Exception(f"gmt file not found: {ontology}")
    else:
        raise Exception(f"Invalid ontology: {ontology}")

    if not isinstance(ontology, list) and ontology.endswith(".gmt"):
        print(f"Running GSEApy with provided gmt file: {ontology}")
        res = run_gseapy(input, ontology, outdir, **kwargs)
        res_merged = res.res2d
        gmt_df = read_gmt(ontology)
        res_merged = res_merged.merge(gmt_df[['Description']], left_index=True, right_index=True, how='left')

    elif ont_id == "GO":
        res_list = []
        for ont in ontology:
            print(f"Running GSEApy with Enrichr library: {ont}")
            res = run_gseapy(input, ont, outdir, **kwargs)
            res_list.append(res.res2d)

        if len(res_list) > 0:
            res_merged = pd.concat(res_list)
            res_merged = res_merged.reset_index(drop=True)

    else:
        print(f"Running GSEApy with Enrichr library: {ontology}")
        res = run_gseapy(input, ontology, outdir, **kwargs)
        res_merged = res.res2d

    convert_gseapy_table(res_merged, ont_id)
    res_merged.to_csv(outfile)


def run_gseapy(input: Union[pd.DataFrame, pd.Series], ontology: str, outdir: str, **kwargs):

    res = gseapy.prerank(rnk=input, 
                          gene_sets=ontology, 
                          outdir=None, 
                          min_size=10,
                          max_size=500,
                          **kwargs)
    res.res2d["Ontology"] = ontology
    return res

def read_gmt(gmt_file):
    gene_sets = []
    with open(gmt_file) as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                description = parts[1]
                genes = parts[2:]
                gene_sets.append({'Description': description, 'genes': genes})
    return pd.DataFrame(gene_sets)

def main() -> None:
    tab = pd.read_csv(input_file, index_col=0)

    if tab.index.name != "SYMBOL":
        gene_table = pd.read_csv(gene_table_file, index_col=0)
        tab[keytype] = tab.index
        tab = tab.merge(gene_table, how='left', on=keytype)
        tab.dropna(axis=0, inplace=True)
        tab.set_index("SYMBOL", inplace=True)

    # https://gseapy.readthedocs.io/en/latest/faq.html#q-why-gene-symbols-in-enrichr-library-are-all-upper-cases-for-mouse-fly-fish-worm
    tab.index = tab.index.str.upper() # Enrichr supports only upper case

    run_gseapy_multi(tab, metric=metric, ontology=ontology, organismKEGG=organismKEGG, outfile=outfile)

if __name__ == "__main__":

    input_file = sys.argv[1]
    keytype = sys.argv[2]
    organismKEGG = sys.argv[3]
    gene_table_file = sys.argv[4]
    metric = sys.argv[5]
    ontology = sys.argv[6] # either "GO", "KEGG", Enrichr library, or path to gmt file
    outfile = sys.argv[7]

    main()

