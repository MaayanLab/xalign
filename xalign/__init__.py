import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import requests, sys
import urllib.request
import os
import xalign.file as filehandler
import sys
import platform
import tarfile

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

def build_index(aligner: str, species: str):
    organisms = retrieve_ensembl_organisms()
    if species in organisms:
        print("Download fastq.gz")
        print(filehandler.get_data_path())
        filehandler.download_file(organisms[species][2], "temp.fastq.gz")

        print("Download gtf")
        filehandler.download_file(organisms[species][3], "temp.gtf")
    else:
        sys.exit(0)
    
    osys = platform.system()

    if aligner == "kallisto":
        url = ""
        if osys == "windows":
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_windows-v0.46.1.zip"
            filepath = filehandler.download_file(url, "kallisto.zip")
            file = tarfile.open(filepath)
            file.extractall('kallisto/kallisto.exe', filehandler.get_data_path())
            file.close()
        elif osys == "linux":
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz"
            filepath = filehandler.download_file(url, "kallisto.tar.gz")
            file = tarfile.open(filepath)
            file.extractall('kallisto/kallisto', filehandler.get_data_path())
            file.close()
        else: #mac
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_mac-v0.46.1.tar.gz"
            filepath = filehandler.download_file(url, "kallisto.tar.gz")
            file = tarfile.open(filepath)
            file.extractall('kallisto/kallisto', filehandler.get_data_path())
            file.close()
        filehandler.download_file()
    elif aligner == "hisat2":
        print("missing")
    elif aligner == "salmon":
        print("missing")