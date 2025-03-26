# scripts/rung_string.py

import os
import sys
import requests
import json
import time
from pathlib import Path

import pandas as pd
from typing import Optional


def check_api_key(key: Optional[str]) -> str:

    if key in ["", "_"] or key is None:
        print("No STRING API key supplied. Enter key or request one by entering y")
        answer = input()
        if answer.lower() == "y":
            print("Requesting API key from STRING")
        elif len(answer) > 1:
            return answer
        else:
            raise Exception("No STRING API key found or requested. Abborting...")
    else:
        return key

    post = "https://version-12-0.string-db.org/api/json/get_api_key"
    response = requests.post(post)
    r = response.content.decode("utf-8")
    key = r.split('api_key": "')[1].split('"')[0]
    note = r.split('note": "')[1].split('"')[0]

    if key is None:
        raise Exception("Error requesting API key from STRING. Abborting...")

    print("API Key:", key)
    print(note)
    print("Please save your key securely. The key is completely anonymous, and if you lose it, we will not be able to "
    "identify you or recover it. Please avoid requesting multiple keys to ensure fair usage of the resources.")

    answer = ""
    while answer.lower() != "y":
        print("Press y to confirm that you have saved your key securely.")
        answer = input()

    return key


def format_string_table(df: pd.DataFrame, library: str) -> pd.DataFrame:
    """
    Format table from STRING databse functional scoring results (proteins with values/ranks)
    Output will look closer to ClusterProfiler table
    """

    df["ONTOLOGY"] = df.index.to_series()

    match library:
        case "GO":
            df = df[df.index.str.startswith("GO ")]
            df = df.replace({"ONTOLOGY": "GO Process"}, "BP")
            df = df.replace({"ONTOLOGY": "GO Function"}, "MF")
            df = df.replace({"ONTOLOGY": "GO Component"}, "CC")
        case "KEGG":
            df = df[df.index.str.startswith("KEGG")]
        case _:
            raise Exception("Unknown library:", library)

    df = df.rename(
        {
            "enrichment score": "enrichmentScore",
            "term description": "Description",
            "term ID": "ID",
            "false discovery rate": "qvalue",
        },
        axis=1
    )
    df = df.rename(
        {
            "enrichmentScore": "enrichmentScore",
            "termDescription": "Description",
            "termID": "ID",
            "falseDiscoveryRate": "qvalue",
        },
        axis=1,
    )
    df.set_index("ID", inplace=True)

    df["pvalue"] = df["qvalue"]  # dubious but STRING doesn't save pvalues...

    # STRING sort values from negative to positive, hence "top" will be downregulated, hence reverse this here
    df["enrichmentScore"] = df["enrichmentScore"] * df["direction"].apply(
        lambda x: -1 if x == "top" else 1 if x == "bottom" else 0
    )

    return df


def prepare_string_input(path: str, metric: str) -> None:
    tab = pd.read_csv(path, index_col=0)
    tab = tab[metric]
    formatted_path = Path(path.replace(".csv", ".string.tsv"))
    tab.to_csv(formatted_path, sep="\t", header=False)


def run_string_enrichment(input_path: str, api_key: str, species: int, fdr: float = 0.05) -> str | None:
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "json"
    method = "valuesranks_enrichment_submit"

    # Construct the request

    request_url = "/".join([string_api_url, output_format, method])

    # Read the input data as string

    identifiers = open(input_path).read()

    # Set parameters

    params = {
        "species": species,  # NCBI/STRING species identifier (e.g., 9606 for human)
        "caller_identity": "www.awesome_app.org",
        "identifiers": identifiers,
        "api_key": api_key,
        "ge_fdr": fdr,
        "ge_enrichment_rank_direction": -1
    }

    # Call STRING

    response = requests.post(request_url, data=params)

    # Read and parse the result

    data = json.loads(response.text)[0]

    if 'status' in data and data['status'] == 'error':
        print("Status:", data['status'])
        print("Message:", data['message'])
        raise Exception("Abborting...")
    else:
        job_id = data["job_id"]
        print(f"Job submitted successfully. Job ID: {job_id}")
        return job_id


def check_job_status(api_key: str, job_id: str) -> dict:

    url = f"https://version-12-0.string-db.org/api/json/valuesranks_enrichment_status?api_key={api_key}&job_id={job_id}"
    response = requests.post(url)
    response_dict = json.loads(response.text)[0]
    return response_dict


def main(input_path: str, key: str, metric: str, outfile: str, fdr: float = 0.05) -> None:

    api_key = check_api_key(key)

    prepare_string_input(input_path, metric)

    input_path = input_path.replace(".csv", ".string.tsv")

    match organism_kegg:
        case "hsa":
            species = 9606
        case "mmu":
            species = 10090
        case _:
            raise Exception(f"Organism not supported: {organism_kegg}")

    outfile_response = outfile.replace(".csv", ".response.json").replace("_PLACEHOLDER_.", "")
    print(outfile_response)
    if os.path.isfile(outfile_response):
        with open(outfile_response) as f:
            response_dict = json.load(f)
        if "job_id" in response_dict:
            job_id = response_dict["job_id"]

            while True:
                print(f"Found existing job_id: {job_id}. Send new STRING Job? [y / n]")
                res = input()
                if res.lower() == "y":
                    job_id = run_string_enrichment(input_path, api_key, species, fdr)
                    break
                elif res.lower() == "n":
                    print("Using existing job id")
                    break

    else:
        job_id = run_string_enrichment(input_path, api_key, species, fdr)

    if job_id is None:
        raise Exception("No job ID found. Abborting...")

    time_waited = 0
    while time_waited < 3600:
        if time_waited % 20 == 0:
            response_dict = check_job_status(api_key, job_id)

            if "status" in response_dict:
                with open(outfile_response, 'w', encoding='utf-8') as f:
                    json.dump(response_dict, f, ensure_ascii=False, indent=4)
                if response_dict["status"] == "success":
                    break
            elif "status" in response_dict and response_dict["status"] == "failed":  # TO DO: check if this is correct
                raise Exception("STRING analysis failed. Abborting...")

        time.sleep(1)
        time_waited += 1
        if response_dict and 'message' in response_dict:
            print(f"\rTime elapsed: {time_waited}s, status: {response_dict['message']}", end='', flush=True)

    else:
        raise Exception("STRING analysis not completed after 3600 seconds. Abborting...")

    download_url = response_dict['download_url']
    df = pd.read_csv(download_url, sep='\t', index_col=0)

    for library in ["KEGG", "GO"]:
        outfile_lib = outfile.replace("_PLACEHOLDER_", library)
        df_lib = format_string_table(df, library=library)
        df_lib.to_csv(outfile_lib)


if __name__ == "__main__":
    input_file = sys.argv[1]
    api_key = sys.argv[2]
    keytype = sys.argv[3]
    organism_kegg = sys.argv[4]
    gene_table_file = sys.argv[5]
    metric = sys.argv[6]
    outfile = sys.argv[7]
    fdr = float(sys.argv[8])

    main(input_file, api_key, metric, outfile, fdr)
