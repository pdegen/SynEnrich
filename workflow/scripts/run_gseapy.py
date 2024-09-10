# scripts/run_gseapy.py

import sys
import numpy as np
import pandas as pd
import gseapy
from typing import Union

def convert_gseapy_table(tab: pd.DataFrame, ont_id: str) -> None:
    '''Convert GSEApy columns to match clusterProfiler'''
    tab.rename({"NOM p-val": "pvalue", "ES": "enrichmentScore", 
                "FDR q-val": "qvalue"}, axis=1, inplace=True) # dubious combining scores?

    if gseapy_use_enrichr or ont_id == "KEGG":
        tab.rename({"Term": "Description"}, axis=1, inplace=True)

    if ont_id == "GO" and gseapy_use_enrichr:
        tab["ID"] = "GO:" + tab["Description"].str.split("\\(GO:").str[1].str[:-1]
    elif ont_id == "GO" and not gseapy_use_enrichr:
        tab.rename({"Term":"ID"}, axis=1, inplace=True)
    elif ont_id == "KEGG":
        tab["ID"] = tab["Description"] # gseapy doesn't save KEGG id, hence use
    else:
        raise Exception("Unknown ontology:", ont_id)
    
    tab.set_index("ID", inplace=True, drop=False)

def run_gseapy_multi(tab: pd.DataFrame, metric: str, outdir: str = "", overwrite: bool = False, **kwargs) -> None:

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

        print("Running:", onts)

        if gseapy_use_enrichr or ont_id != "GO":
            res_list = []
            for ont in onts:
                print("Running GSEApy with Enrichr libraries:", ont)
                res = run_gseapy(input, ont, outdir, **kwargs)
                res_list.append(res.res2d)

            res_merged = pd.concat(res_list)
            res_merged = res_merged.reset_index(drop=True)

        else:
            print("Running GSEApy with provided gmt files")
            if organismKEGG == "hsa":
                gmt_file = "resources/go_terms.org.Hs.eg.db.gmt"
            elif organismKEGG == "mmu":
                gmt_file = "resources/go_terms.org.Mm.eg.db.gmt"

            res = run_gseapy(input, gmt_file, outdir, **kwargs)
            res_merged = res.res2d
            gmt_df = read_gmt(gmt_file)
            res_merged = res_merged.merge(gmt_df[['Description']], left_index=True, right_index=True, how='left')

        convert_gseapy_table(res_merged, ont_id)
        res_merged.to_csv(outfile_go if ont_id == "GO" else outfile_kegg)


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

    run_gseapy_multi(tab, metric)

if __name__ == "__main__":

    input_file = sys.argv[1]
    keytype = sys.argv[2]
    organismKEGG = sys.argv[3]
    gene_table_file = sys.argv[4]
    metric = sys.argv[5]
    outfile_go = sys.argv[6]
    outfile_kegg = sys.argv[7]
    gseapy_use_enrichr = sys.argv[8]
    gseapy_use_enrichr = gseapy_use_enrichr == "True"

    ontologies = ["GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"]

    if organismKEGG == "hsa":
        ontologies += ["KEGG_2021_Human"]
    elif organismKEGG == "mmu":
        ontologies += ["KEGG_2019_Mouse"]
    else:
        raise Exception(f"Organism not yet supported: {organismKEGG}")
    main()

