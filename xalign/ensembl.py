import numpy as np
import pandas as pd
import requests, sys
import os
import xalign.file as filehandler
import sys
import json
import mygene
import requests, sys
import json

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
        gtf_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/gtf/"+name+"/"+name.capitalize()+"."+assembly+"."+str(release)+".gtf.gz"
        organisms[name] = [name, disp, cdna_url, gtf_url]
        
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

def map_ensembl_genesymbol(ids):
    mg = mygene.MyGeneInfo()
    return mg.querymany(ids, scopes='ensembl.transcript')

def agg_gene_counts(transcript_counts, species, identifier="symbol"):
    
    transcript_counts.index = transcript_counts.iloc[:, 0].str.replace("\.[0-9]", "", regex=True)
    print(transcript_counts)
    
    if not os.path.exists(species+"_ensembl_ids.json"):
        ids = list(transcript_counts.index)
        mg = mygene.MyGeneInfo()
        id_query = mg.querymany(ids, scopes='ensembl.transcript', fields=["ensembl", "symbol", "entrezgene", "name"])

        jd = json.dumps(id_query)
        f = open(species+"_ensembl_ids.json","w")
        f.write(jd)
        f.close()
    else:
        f = open(species+"_ensembl_ids.json","r")
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
            ensemblid = q["ensembl"]["gene"]
        ginfo.append([q["query"], symbol, ensemblid, entrezgene, name])

    gene_map = pd.DataFrame(ginfo)
    gene_map.index = gene_map.iloc[:,0]
    
    tc = transcript_counts.join(gene_map, how="inner")
    tc.columns = ["transcript", "counts", "tpm", "transcript2", "symbol", "ensembl_id", "entrezgene_id", "name"]
    
    tc = tc.groupby([identifier], as_index=False)['est_counts'].agg('sum')
    tc.iloc[:,1] = tc.iloc[:,1].astype("int")
    
    return tc[tc[identifier] != ""]