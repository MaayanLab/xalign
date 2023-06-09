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
import biomart

import xalign.file as filehandler

def string_vector_to_int_vector(string_vector):
    unique_elements = sorted(list(set(string_vector)))
    int_vector = [unique_elements.index(elem) for elem in string_vector]
    return int_vector

def download_and_process_sample(sample, sra_files, ensembl_idx):
    gene_counts = []
    for fname in sra_files:
        tmp_file = ""
        try:
            kallisto_result = pd.read_csv(tmp_file.name, sep="\t")
            values = kallisto_result.iloc[:, 1].to_numpy(dtype=float)
            gene_counts.append(np.bincount(ensembl_idx, weights=values))
        except Exception:
            gene_counts.append([0] * len(set(ensembl_idx)))
        finally:
            os.remove(tmp_file.name)
    max_value_uint32 = np.iinfo(np.uint32).max
    return np.clip(np.array(gene_counts).sum(axis=0), 0, max_value_uint32).astype(np.uint32)

def get_ensembl_mappings(species):
    string = species.split('_')
    kk = ""
    for s in string[:-1]:
        kk = kk+s[0]
    kk+string[-1]
    # Set up connection to server
    #server = biomart.BiomartServer('http://useast.ensembl.org/biomart')
    #version 2 (107)
    server = biomart.BiomartServer('http://jul2022.archive.ensembl.org/biomart')
    if species == "mouse":
        mart = server.datasets['mmusculus_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'mgi_symbol', 'ensembl_gene_id', 'gene_biotype']
    else:
        mart = server.datasets['hsapiens_gene_ensembl']
        attributes = ['ensembl_transcript_id', 'hgnc_symbol', 'ensembl_gene_id', 'gene_biotype']                                                     
    # Get the mapping between the attributes                                    
    response = mart.search({'attributes': attributes})
    data = response.raw.data.decode('ascii')
    ensembl_ids = []
    # Store the data in a dict                                                  
    for line in data.splitlines():                                              
        line = line.split('\t')                                
        ensembl_ids.append(line)
    gene_map = pd.DataFrame(ensembl_ids)
    gene_map.index = gene_map.iloc[:,0]
    nn = np.where(gene_map.iloc[:,1] == "")[0]
    gene_map.iloc[nn, 1] = gene_map.iloc[nn, 2]
    gene_map.columns = ["ensembl_transcript", "symbol", "ensembl_gene", "biotype"]
    gene_map = gene_map[~gene_map.index.duplicated(keep='first')]
    return gene_map

def retrieve_ensembl_organisms(release=None):
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
        if not release:
            release = sp["release"]
        name = sp["name"]
        disp = sp["display_name"]
        assembly = sp["assembly"]
        cdna_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/fasta/"+name+"/cdna/"+name.capitalize()+"."+assembly+".cdna.all.fa.gz"
        ncdna_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/fasta/"+name+"/ncrna/"+name.capitalize()+"."+assembly+".ncrna.fa.gz"
        gtf_url = "http://ftp.ensembl.org/pub/release-"+str(release)+"/gtf/"+name+"/"+name.capitalize()+"."+assembly+"."+str(release)+".gtf.gz"
        organisms[name] = [name, disp, cdna_url, gtf_url, ncdna_url, release]
        
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
