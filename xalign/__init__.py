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

def build_index(aligner: str, species: str, overwrite=False):
    organisms = retrieve_ensembl_organisms()
    if species in organisms:
        print("Download fastq.gz")
        print(filehandler.get_data_path())
        filehandler.download_file(organisms[species][2], "temp.fastq.gz")
        #print("Download gtf")
        #filehandler.download_file(organisms[species][3], "temp.gtf")
    else:
        sys.exit(0)
    osys = platform.system().lower()
    download_aligner(aligner, osys)
    os.makedirs(filehandler.get_data_path()+"/index", exist_ok=True)
    if aligner == "kallisto":
        if (not os.path.exists(filehandler.get_data_path()+"index/"+species+".idx")) or overwrite:
            print("Build kallisto index for "+species)
            os.system(filehandler.get_data_path()+"/kallisto/kallisto index -i "+filehandler.get_data_path()+"index/"+species+".idx "+filehandler.get_data_path()+"temp.fastq.gz")
    else:
        print("Index already exists. Use overwrite to rebuild.")

def download_aligner(aligner, osys):
    
    print(osys)

    if aligner == "salmon":
        if osys == "linux":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extract('salmon-1.5.2_linux_x86_64/bin/salmon', filehandler.get_data_path())
            file.close()
        else:
            print("Salmon not supported by this package for this operating system.")

    if aligner == "kallisto":
        if osys == "windows":
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_windows-v0.46.1.zip"
            filepath = filehandler.download_file(url, "kallisto.zip")
            file = tarfile.open(filepath)
            file.extract('kallisto/kallisto.exe', filehandler.get_data_path())
            file.close()
        elif osys == "linux":
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz"
            filepath = filehandler.download_file(url, "kallisto.tar.gz")
            file = tarfile.open(filepath)
            file.extract('kallisto/kallisto', filehandler.get_data_path())
            file.close()
        else: #mac
            url = "https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_mac-v0.46.1.tar.gz"
            filepath = filehandler.download_file(url, "kallisto.tar.gz")
            file = tarfile.open(filepath)
            file.extract('kallisto/kallisto', filehandler.get_data_path())
            file.close()
    elif aligner == "hisat2":
        print("missing")
    elif aligner == "salmon":
        print("missing")
