# scripts/run_gseapy.py

import os
import sys
import yaml
import numpy as np
import pandas as pd
import gseapy
from typing import Union

def convert_gseapy_table(tab: pd.DataFrame, ont_id: str) -> None:
    '''Convert GSEApy columns to match clusterProfiler'''
    tab.rename({"NOM p-val": "pvalue", "ES": "enrichmentScore", 
                "Term": "Description", "FDR q-val": "qvalue"}, axis=1, inplace=True) # dubious combining scores?

    if ont_id == "GO":
        tab["ID"] = "GO:" + tab["Description"].str.split("\\(GO:").str[1].str[:-1]
    elif ont_id == "KEGG":
        tab["ID"] = tab["Description"] # gseapy doesn't save KEGG id, hence use description
    else:
        raise Exception("Unknown ontology:", ont_id)
    
    tab.set_index("ID", inplace=True, drop=False)

def run_gseapy_prerank_multi(tab: pd.DataFrame, 
                     metric: str, 
                     outdir: str = None, 
                     overwrite: bool = False, 
                     **kwargs) -> None:

    if metric not in tab.columns:
        if metric == "neg_signed_logpval":
            tab["neg_signed_logpval"] = -np.sign(tab["logFC"]) * np.log10(tab["PValue"])
        else:
            raise Exception(f"Metric not in input table: {metric}")
        
    input = tab[metric].sort_values(ascending=False)
    input.index.name = None # remove header

    go_onts = [o for o in ontologies if "GO" in o and "KEGG" not in o]
    kegg_onts = [o for o in ontologies if "KEGG" in o]

    for onts in [go_onts, kegg_onts]:
        ont_id = "GO" if onts == go_onts else "KEGG"

        res_list = []
        # TO DO: check if looping increases runtime compared to passing list of onts; currently, we loop to do res.res2d["Ontology"] = ontology
        for ont in onts:
            print("Running prerank ont:", ont)
            res = run_gseapy_prerank(input, ont, outdir, **kwargs)
            res_list.append(res.res2d)

        # TO DO: think about whether FDR correction should be re-applied on merged list
        res_merged = pd.concat(res_list)
        res_merged = res_merged.reset_index(drop=True)

        convert_gseapy_table(res_merged, ont_id)
        res_merged.to_csv(outfile_go if ont_id == "GO" else outfile_kegg)


def run_gseapy_prerank(input: Union[pd.DataFrame, pd.Series], 
               ontology: str, 
               outdir: str = None,
               **kwargs) -> gseapy.Prerank:
    
    res = gseapy.prerank(rnk=input, 
                          gene_sets=ontology, 
                          outdir=outdir, 
                          min_size=10,
                          max_size=500,
                          **kwargs)
    res.res2d["Ontology"] = ontology
    return res

def run_gseapy_multi(input: pd.DataFrame, 
                     metric: str,
                     pheno_pos: str,
                     pheno_neg: str,
                     **kwargs) -> None:
        

    go_onts = [o for o in ontologies if "GO" in o and "KEGG" not in o]
    kegg_onts = [o for o in ontologies if "KEGG" in o]

    for onts in [go_onts, kegg_onts]:
        ont_id = "GO" if onts == go_onts else "KEGG"

        res_list = []
        # TO DO: check if looping increases runtime compared to passing list of onts; currently, we loop to do res.res2d["Ontology"] = ontology
        for ont in onts:
            print("Running standard ont:", ont)
            res = run_gseapy(input, ont, metric, pheno_pos, pheno_neg, **kwargs)
            res_list.append(res.res2d)

        # TO DO: think about whether FDR correction should be re-applied on merged list
        res_merged = pd.concat(res_list)
        res_merged = res_merged.reset_index(drop=True)

        convert_gseapy_table(res_merged, ont_id)
        res_merged.to_csv(outfile_go if ont_id == "GO" else outfile_kegg)

def run_gseapy(input: pd.DataFrame,
               ontology: str,
               metric: str,
               pheno_pos: str,
               pheno_neg: str,
               **kwargs) -> gseapy.GSEA:
    
    gs = gseapy.GSEA(data=input,
            gene_sets=ontology,
            classes = class_vector, # cls=class_vector
            permutation_type='gene_set' if len(input.columns) < 15 else "phenotype",
            permutation_num=1000, # reduce number to speed up test
            outdir=None,
            min_size = 10,
            max_size = 500,
            method=metric)
    gs.pheno_pos = pheno_pos
    gs.pheno_neg = pheno_neg
    gs.run()
    return gs

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
    
    if is_prerank:
        run_gseapy_prerank_multi(tab, metric)
    else:
        tab = tab.iloc[:,:len(class_vector)] # remove ENTREZID
        run_gseapy_multi(tab, metric, pheno_pos, pheno_neg)


if __name__ == "__main__":


    config_file_path = os.path.join("config", "config.yaml")
    with open(config_file_path, 'r') as stream:
        config = yaml.safe_load(stream)

    class_vector = config.get('class_vector', [])
    metrics_prerank = config.get('metrics_prerank', [])
    organismKEGG = config.get('organismKEGG')
    keytype = config.get('keytype')
    pheno_neg = config.get('pheno_pos')
    pheno_pos = config.get('pheno_neg')

    input_file: str = sys.argv[1]
    gene_table_file: str = sys.argv[2]
    metric: str = sys.argv[3]
    outfile_go: str = sys.argv[4]
    outfile_kegg: str = sys.argv[5]

    is_prerank = metric in metrics_prerank

    ontologies = ["GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"]

    if organismKEGG == "hsa":
        ontologies += ["KEGG_2021_Human"]
    elif organismKEGG == "mmu":
        ontologies += ["KEGG_2019_Mouse"]
    else:
        raise Exception(f"Organism not yet supported: {organismKEGG}")
    main()

