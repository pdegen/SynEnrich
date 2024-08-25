import pickle
import pandas as pd

def pickler(contents, filepath):
    """Store arbitrary Python object in filepath"""
    with open(filepath, "wb") as fp:
        pickle.dump(contents, fp)

def format_string_table(df: pd.DataFrame, library: str) -> pd.DataFrame:
    """
    Format table from STRING databse functional scoring results (proteins with values/ranks)
    Output will look closer to ClusterProfiler table
    """

    df.loc[:,"ONTOLOGY"] = df.index

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