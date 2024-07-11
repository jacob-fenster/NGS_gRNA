import requests
import time
import xml.etree.ElementTree as ET
import pdb
import sys
import ast

def request_data(url, max_retries=100, delay=0.11):
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response
            else:
                retries += 1
                time.sleep(delay)
        except:
            retries += 1
            time.sleep(delay)
    print(f"Max retries of {max_retries} reached with url:\n{url}\nExiting")
    sys.exit()

def check_200_request_error(url, response, max_retries=200):
        retries = 0
        while retries < max_retries:
            try:
                root = ET.fromstring(response.text)
                while root.find('./ERROR') != None and retries < max_retries:
                    response = request_data(url)
                    try:
                        root = ET.fromstring(response.text)
                    except: #no automated logic here yet
                        retries += 1
                if retries < max_retries:
                    return response
                else:
                    print(f"Max retries of {max_retries} reached with url below in check_200_request_error:\n{url}\nExiting")
                    sys.exit()
            except: #no automated logic here yet
                retries += 1
                response = request_data(url)
    
def esearch(query, db, api_key):
    base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    url = f"{base}?db={db}&term={query}&api_key={api_key}&usehistory=y"
    request = request_data(url)
    request = check_200_request_error(url, request)
    try: #try and extract the Count information
        root = ET.fromstring(request.text)
        max_count, QueryKey, WebEnv = int(root.find('./Count').text), root.find('./QueryKey').text, root.find('./WebEnv').text
    except: #this is here to troublshoot. Nothing automated to handle exceptions yet
        pdb.set_trace()
        x = 1
    print(f"There are {str(max_count)} results retured for query:\n {query.replace('+', ' ')}")
    return QueryKey, WebEnv, max_count

def efetch_aaCDS_text(QueryKey, WebEnv, max_count, api_key, output_basename):
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    db, rettype, retmode, retmax = 'nuccore', 'fasta_cds_aa', 'text', 200
    #clear files with at start
    with open(f"{output_basename}.fasta", 'w') as f:
        pass
    with open(f"{output_basename}_fetch_summary.tsv", 'w') as f:
        f.write(f"#QueryKey={QueryKey}, WebEnv={WebEnv}, IDs={max_count}\n")
    with open(f"{output_basename}_failed_fetch.txt", "w") as f:
        pass
    for retstart in range(0, max_count, retmax):
        url=f"{base}?db={db}&query_key={QueryKey}&WebEnv={WebEnv}&retstart={retstart}&rettype={rettype}&rettmode={retmode}&retmax={retmax}&api_key={api_key}"
        fetch = request_data(url)
        if retstart % 100000 == 0:
            print(f"Efetch {int((retstart/max_count)*100)}%. query {retstart} of {max_count}")
        if fetch.headers['Content-Type'] == 'text/plain': #see if correct text output format and contains information
            handle = open(f"{output_basename}.fasta", 'a')
            summary = open(f"{output_basename}_fetch_summary.tsv", 'a')
            #filter out blank lines returning a list of stripped
            fetch_lines = fetch.text.split("\n")
            filtered_fetch = [line for line in fetch_lines if line.strip()]
            if len(filtered_fetch) == 0:
                first_line = 'No Data'
            else:
                first_line = filtered_fetch[0]
            #can do stats here if needed
            summary.write(f"retstart:{retstart}\tfirst_line:{first_line}\tPrefiltered_lines:{len(fetch_lines)}\tFiltered_lines:{len(filtered_fetch)}\n")
            for line in filtered_fetch:
                handle.write(f"{line}\n")
            handle.close()
            summary.close()
        else:#write error file if didn't get clusters
            print(f"Fetch failed with status_code={str(fetch.status_code)}\tretstart={retstart}\tretmax={retmax}\tquery_key={QueryKey}\tWebEnv={WebEnv}\nText={fetch.text}\n")
            f = open(f"{output_basename}_failed_fetch.txt", "a") 
            f.write(f"status_code={str(fetch.status_code)}\tretstart={retstart}\tretmax={retmax}\tquery_key={QueryKey}\tWebEnv={WebEnv}\nText={fetch.text}\n")
            f.close()
    print("Efetch succesfully finished")
api_key = "46f09d13b8f485ab42dabf140b3a09dd9c09"

query_all_NCBI = 'Viruses[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])'
query_all_virus = 'Viruses[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[PROP]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[FILTER]+OR+complete+genome[TITL])'
query_ASFV = 'Viruses[Organism]+AND+African+swine+fever+virus[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[PROP]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[FILTER]+OR+complete+genome[TITL])'
output_basename_all_NCBI = '/Users/jacobfenster/Local_Documents/Databases/proteins_of_all_viral_genomes/NCBI_all_viral_genomes_aaCDS_v2'
output_basename_all = '/Users/jacobfenster/Local_Documents/Databases/proteins_of_all_viral_genomes/CompleteGenome_virus_CDS_NCBI_v2'
output_basename_ASFV = '/Users/jacobfenster/Local_Documents/Databases/proteins_of_all_viral_genomes/ASFV_all_NCBI'
db = 'nuccore'

QueryKey, WebEnv, max_count = esearch(query_all_NCBI, db, api_key)
efetch_aaCDS_text(QueryKey, WebEnv, max_count, api_key, output_basename_all_NCBI)

