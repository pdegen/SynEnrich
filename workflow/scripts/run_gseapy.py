# scripts/run_gseapy.py

import sys
import pandas as pd
import gseapy

def convert_gseapy_table(tab: pd.DataFrame, ont_id: str) -> None:
    '''Convert GSEApy columns to match clusterProfiler'''
    tab.rename({"NOM p-val": "pvalue", "ES": "enrichmentScore", 
                "Term": "Description", "FDR q-val": "qvalue"}, axis=1, inplace=True) # dubious combining scores?

    if ont_id == "GO":
        tab["ID"] = "GO:" + tab["Description"].str.split("\\(GO:").str[1].str[:-1]
    elif ont_id == "KEGG":
        tab["ID"] = tab["Description"] # gseapy doesn't save KEGG id, hence use
    else:
        raise Exception("Unknown ontology:", ont_id)
    
    tab.set_index("ID", inplace=True, drop=False)

def run_gseapy_multi(tab, metric, outdir=None, overwrite=False, **kwargs):

    input = tab[metric].sort_values(ascending=False)
    input.index.name = None # remove header

    go_onts = [o for o in ontologies if "GO" in o and "KEGG" not in o]
    kegg_onts = [o for o in ontologies if "KEGG" in o]

    for onts in [go_onts, kegg_onts]:
        ont_id = "GO" if onts == go_onts else "KEGG"

        res_list = []
        print("Running:", onts)
        for ont in onts:
            print("Sub-Running:", ont)
            res = run_gseapy(input, ont, outdir, **kwargs)
            res_list.append(res.res2d)

        res_merged = pd.concat(res_list)
        res_merged = res_merged.reset_index(drop=True)

        convert_gseapy_table(res_merged, ont_id)
        print(res_merged.head())
        res_merged.to_csv(outfile_go if ont_id == "GO" else outfile_kegg)


def run_gseapy(input, ontology, outdir=None, **kwargs):
    res = gseapy.prerank(rnk=input, 
                          gene_sets=ontology, 
                          outdir=None, 
                          min_size=10,
                          max_size=500,
                          **kwargs)
    res.res2d["Ontology"] = ontology
    return res

def main():
    tab = pd.read_csv(input_file, index_col=0)

    if tab.index.name != "SYMBOL":
        gene_table = pd.read_csv(gene_table_file, index_col=0)
        tab[keytype] = tab.index
        tab = tab.merge(gene_table, how='left', on=keytype)
        tab.dropna(axis=0, inplace=True)
        tab.set_index("SYMBOL", inplace=True)
    else:
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

    ontologies = ["GO_Biological_Process_2023","GO_Cellular_Component_2023","GO_Molecular_Function_2023"]

    if organismKEGG == "hsa":
        ontologies += ["KEGG_2021_Human"]
    elif organismKEGG == "mmu":
        ontologies += ["KEGG_2019_Mouse"]
    else:
        raise Exception(f"Organism not yet supported: {organismKEGG}")
    main()

