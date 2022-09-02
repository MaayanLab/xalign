import numpy as np
import pandas as pd
import requests, sys
import os
import sys
import json
import mygene
import requests, sys
import json
import multiprocessing
from tqdm import tqdm
from itertools import chain


import xalign.file as filehandler

def retrieve_ensembl_organisms():
    server = "http://rest.ensembl.org"
    ext = "/info/species?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    species = decoded["species"]
    organisms = {}
    
    for sp in species:
        release = sp["release"]
        name = sp["name"]
        disp = sp["display_name"]
        assembly = sp["assembly"]
        cdna_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/fasta/"+name+"/cdna/"+name.capitalize()+"."+assembly+".cdna.all.fa.gz"
        ncdna_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/fasta/"+name+"/ncrna/"+name.capitalize()+"."+assembly+".ncrna.fa.gz"
        gtf_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/gtf/"+name+"/"+name.capitalize()+"."+assembly+"."+str(release)+".gtf.gz"
        organisms[name] = [name, disp, cdna_url, gtf_url, ncdna_url]
        
    return organisms

def organism_display_to_name(display_name):
    server = "http://rest.ensembl.org"
    ext = "/info/species?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    species = decoded["species"]
    
    for sp in species:
        if display_name == sp["display_name"]:
            return sp["name"]

    return "missing"

def chunk(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def retrieve_ensemble_ids(ids):
    
    chunked_ids = chunk(ids, 1000)
    transcript_info = {}

    server = "https://rest.ensembl.org"
    ext = "/lookup/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    counter = 1
    for cids in chunked_ids:
        r = requests.post(server+ext, headers=headers, data=json.dumps({ "ids" : cids }))
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        print(counter)
        counter = counter + 1
        transcript_info.update(r.json())

def map_transcript(ids):
    mg = mygene.MyGeneInfo()
    return mg.querymany(ids, scopes='ensembl.transcript', fields=["ensembl", "symbol", "entrezgene", "name"], verbose=False)

def chunk(l, n):
	return [l[i:i+n] for i in range(0, len(l), n)]

def agg_gene_counts(transcript_counts, species, identifier="symbol", overwrite=False):
    
    transcript_counts.index = transcript_counts.iloc[:, 0].str.replace("\.[0-9]", "", regex=True)
    
    if not os.path.exists(filehandler.get_data_path()+species+"_ensembl_ids.json") or overwrite:
        ids = list(transcript_counts.index)
        cids = chunk(ids, 200)
        with multiprocessing.Pool(8) as pool:
	        res = list(tqdm(pool.imap(map_transcript, cids), desc="Mapping transcripts", total=len(cids)))
        id_query = list(chain.from_iterable(res))
        jd = json.dumps(id_query)
        f = open(filehandler.get_data_path()+species+"_ensembl_ids.json","w")
        f.write(jd)
        f.close()
    else:
        f = open(filehandler.get_data_path()+species+"_ensembl_ids.json","r")
        id_query = json.load(f)
        f.close()
    
    ginfo = []

    for q in id_query:
        symbol = ""
        entrezgene = ""
        ensemblid = ""
        name = ""
        if "symbol" in q.keys():
            symbol = q["symbol"]
        if "entrezgene" in q.keys():
            entrezgene = q["entrezgene"]
        if "name" in q.keys():
            name = q["name"]
        if "ensembl" in q.keys():
            if isinstance(q["ensembl"], list):
                for x in q["ensembl"]:
                    if x["transcript"] == q["query"]:
                        ensemblid = x["gene"]
            else:         
                ensemblid = q["ensembl"]["gene"]
        ginfo.append([q["query"], symbol, ensemblid, entrezgene, name])

    gene_map = pd.DataFrame(ginfo)
    gene_map.index = gene_map.iloc[:,0]
    
    tc = transcript_counts.join(gene_map, how="inner")
    tc.columns = ["transcript", "counts", "tpm", "transcript2", "symbol", "ensembl_id", "entrezgene_id", "name"]
    
    tc = tc.groupby([identifier], as_index=False)['counts'].agg('sum')
    tc.iloc[:,1] = tc.iloc[:,1].astype("int")
    
    return tc[tc[identifier] != ""]
