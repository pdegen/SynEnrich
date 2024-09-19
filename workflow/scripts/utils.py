import os
import pickle
import pandas as pd
import yaml
from typing import Dict, Any

try:
    import mygene
except ModuleNotFoundError:
    print("mygene not found")
    pass

def pickler(contents, filepath):
    """Store arbitrary Python object in filepath"""
    with open(filepath, "wb") as fp:
        pickle.dump(contents, fp)

def format_string_table(df: pd.DataFrame, library: str) -> pd.DataFrame:
    """
    Format table from STRING databse functional scoring results (proteins with values/ranks)
    Output will look closer to ClusterProfiler table
    """

    df["ONTOLOGY"] = df.index.to_series()

    match library:
        case "GO":
            df = df[df.index.str.startswith("GO ")]
            df.replace({"ONTOLOGY" : "GO Process"}, "BP", inplace=True)
            df.replace({"ONTOLOGY" : "GO Function"}, "MF", inplace=True)
            df.replace({"ONTOLOGY" : "GO Component"}, "CC", inplace=True)
        case "KEGG":
            df = df[df.index.str.startswith("KEGG")]
        case _:
            raise Exception("Unknown library:", library)
        
    df.rename({"enrichment score": "enrichmentScore",
               "term description": "Description",
               "term ID": "ID",
               "false discovery rate": "qvalue"},
              axis=1, inplace=True)
    df.set_index("ID", inplace=True)

    df["pvalue"] = df["qvalue"] # dubious but STRING doesn't save pvalues...

    # STRING sort values from negative to positive, hence "top" will be downregulated, hence reverse this here
    df["enrichmentScore"] = df["enrichmentScore"] * df["direction"].apply(lambda x: -1 if x == "top" else 1 if x == "bottom" else 0)

    return df

def load_config(file_path: str) -> Dict[str, Any]:
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def ensp_to_gene_symbol(ensp_ids, species):
    mg = mygene.MyGeneInfo()
    print("quering for", species)
    # Query mygene.info with the list of ENSP IDs
    result = mg.querymany(ensp_ids, scopes='ensembl.protein', fields='symbol', species=species)

    # Extract gene symbols
    symbols = {item['query']: item.get('symbol', 'N/A') for item in result}
    
    return symbols

def create_gmt(df, output_file):
    # Group by term and description (or category), and collect the ensemblprot IDs into lists
    grouped = df.groupby(['term', 'category', 'description'])['#string_protein_id'].apply(list).reset_index()
    
    with open(output_file, 'w') as f:
        for _, row in grouped.iterrows():
            # Write each line in the required GMT format
            term = row['term']
            cat = row['category']
            description = row['description']
            ensemblprot_list = row['#string_protein_id']
            # Format the row: term, cat, description, and list of ensemblprots
            f.write(f"{term}\t{cat}\t{description}\t" + "\t".join(ensemblprot_list) + "\n")

def create_string_gmt(infile, outfile, orgid, species=""):
    """
    Convert STRING db enrichment file to gmt file
    """
    if species == "": species = orgid

    tab = pd.read_csv(infile, sep="\t")
    tab = tab[tab["term"].str.startswith("GO:")]
    
    prot2symbol_file = f"../../resources/Ontologies/prot2symbol.{species}.csv"

    if os.path.isfile(prot2symbol_file):
        gene_symbols = pd.read_csv(prot2symbol_file, index_col=0)
        gene_symbols = gene_symbols["SYMBOL"].to_dict()
    else:
        print("Retrieving gene symbols...")
        protids = list(set(tab["#string_protein_id"].str.split(f"{orgid}.").str[1]))
        gene_symbols = ensp_to_gene_symbol(protids, species=species)
        df = pd.DataFrame(gene_symbols.values(), index=gene_symbols.keys(), columns=["SYMBOL"])
        df.to_csv(prot2symbol_file)

    tab["#string_protein_id"] = tab["#string_protein_id"].str.replace(f"{orgid}.","")
    tab["#string_protein_id"] = tab["#string_protein_id"].map(gene_symbols).fillna(tab["#string_protein_id"])
    tab["#string_protein_id"] = tab["#string_protein_id"].str.upper()

    cat_dict = {"Biological Process (Gene Ontology)": "BP",
                "Cellular Component (Gene Ontology)": "CC",
                "Molecular Function (Gene Ontology)": "MF"}

    tab["category"] = tab["category"].replace(cat_dict)

    create_gmt(tab, outfile)

def read_enrichr(gmt_file):
    gene_sets = []
    with open(gmt_file) as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                description = parts[0]
                genes = parts[2:] # 2nd col is empty?
                gene_sets.append({'Description': description, 'Genes': genes})
    df = pd.DataFrame(gene_sets)
    df.index = "GO:" + df["Description"].str.split("\(GO:").str[1].str[:-1]
    df.index.name = "ID"
    df["Description"] = df["Description"].str.split("\(GO:").str[0]
    return df