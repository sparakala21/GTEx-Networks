import json
import os
import requests
import mygene
import time

def get_gene_disease_info(HGNC, mg, page=0):
    result = mg.query(HGNC, scopes='symbol', fields='entrezgene', species='human')
    params = {}
    # ...retrieve disease associated to gene with NCBI ID equal to 351
    params["gene_ncbi_id"] = -1
    for i in result['hits']:
        if 'entrezgene' in i:
            params["gene_ncbi_id"] = i['entrezgene']
            break
    if params["gene_ncbi_id"] == -1:
        print(f"Could not find NCBI ID for {HGNC}")
        return None
    # ...retrieve the first page of results (page number 0) 
    params["page_number"] = page

    HTTPheadersDict = {}
    # Set the 'Authorization' HTTP header equal to API_KEY (your API key) 
    HTTPheadersDict['Authorization'] = "aa9215e1-382e-40c4-9880-1bf5e9f512d9"
    # Set the 'accept' HTTP header to specify the response format: one of 'application/json', 'application/xml', 'application/csv' 
    HTTPheadersDict['accept'] = 'application/json'

    response = requests.get("https://api.disgenet.com/api/v1/gda/summary",\
                        params=params, headers=HTTPheadersDict, verify=False)
    if not response.ok:
        if response.status_code == 429:
            while response.ok is False:
                print("You have reached a query limit for your user. Please wait {} seconds until next query".format(\
                response.headers['x-rate-limit-retry-after-seconds']))
                time.sleep(int(response.headers['x-rate-limit-retry-after-seconds']))
                print("Your rate limit is now restored")

                # Repeat your query
                response = requests.get("https://api.disgenet.com/api/v1/gda/summary",\
                                        params=params, headers=HTTPheadersDict, verify=False)
                if response.ok is True:
                    break
                else:
                    continue
    response_parsed = json.loads(response.text)
    if "error" in response_parsed:
        print(f"Error for {HGNC}: {response_parsed['error']}")
        return None
    if "status" not in response_parsed:
        print(f"Unexpected response for {HGNC}: {response_parsed}")
        return None
        print('STATUS: {}'.format(response_parsed["status"]))
        print('TOTAL NUMBER OF RESULTS: {}'.format(response_parsed["paging"]["totalElements"]))
    return response_parsed
if __name__ == "__main__":
    mg = mygene.MyGeneInfo()
    genes = json.load(open("coding_genes.json", "r"))
    num_genes = len(genes)
    i=0
    for gene in genes:
        i+=1
        print(f"Gathering disease {i}/{num_genes}info for {gene}")
        response_parsed = get_gene_disease_info(gene, mg)
        if response_parsed is None:
            continue
        json.dump(response_parsed, open(os.path.join("data", f"{gene}_disease_info.json"), "w"), indent=4)
        genes.remove(gene)
        with open("coding_genes.json", "w") as f:
            json.dump(genes, f, indent=4)


