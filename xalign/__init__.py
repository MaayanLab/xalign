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
import subprocess

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

def build_index(aligner: str, species: str, overwrite=False, verbose=False):
    organisms = retrieve_ensembl_organisms()
    if species in organisms:
        print("Download fastq.gz")
        print(filehandler.get_data_path())
        filehandler.download_file(organisms[species][2], species+".fastq.gz", overwrite=overwrite, verbose=False)
        #print("Download gtf")
        #filehandler.download_file(organisms[species][3], species+".gtf")
    else:
        sys.exit(0)
    osys = platform.system().lower()
    download_aligner(aligner, osys)
    os.makedirs(filehandler.get_data_path()+"/index", exist_ok=True)
    if aligner == "kallisto":
        if (not os.path.exists(filehandler.get_data_path()+"index/kallisto_"+species+".idx")) or overwrite:
            print("Build kallisto index for "+species)
            os.system(filehandler.get_data_path()+"kallisto/kallisto index -i "+filehandler.get_data_path()+"index/kallisto_"+species+".idx "+filehandler.get_data_path()+species+".fastq.gz")
        else:
            print("Index already exists. Use overwrite to rebuild.")
    elif aligner == "salmon":
        if (not os.path.exists(filehandler.get_data_path()+"index/salmon_"+species)) or overwrite:
            print("Build salmon index for "+species)
            print(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin index -i "+filehandler.get_data_path()+"index/salmon_"+species+".idx -t "+filehandler.get_data_path()+species+".fastq.gz")
            os.system(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon index -i "+filehandler.get_data_path()+"index/salmon_"+species+" -t "+filehandler.get_data_path()+species+".fastq.gz")
        else:
            print("Index already exists. Use overwrite to rebuild.")

def download_aligner(aligner, osys):
    
    print(osys)

    if aligner == "salmon":
        if osys == "linux":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extractall(filehandler.get_data_path())
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


def align_fastq(aligner, species, fastq, t=1, overwrite=False, verbose=False):
    if isinstance(fastq, str):
        fastq = [fastq]
    build_index(aligner, species, overwrite=overwrite, verbose=False)
    if aligner == "kallisto":
        if len(fastq) == 1:
            print("Align with kallisto (single strand).")
            res = subprocess.Popen(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto --single -l 200 -s 20 "+fastq[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        else:
            print("Align with kallisto (paired).")
            res = subprocess.Popen(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto "+fastq[0]+" "+fastq[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if res.wait() != 0:
            output, error = res.communicate()
            print(error)
            print(output)
    elif aligner == "salmon":
        if len(fastq) == 1:
            print("Align with salmon (single).")
            res = subprocess.Popen(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/salmon_"+species+" -l A -r "+fastq[0]+" -p "+str(t)+" --validateMappings -o "+filehandler.get_data_path()+"outsalmon", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if res.wait() != 0:
                output, error = res.communicate()
                print(error)
                print(output)
        else:
            print("Align with salmon (single).")
            res = subprocess.Popen(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/salmon_"+species+" -l A -1 "+fastq[0]+" -2 "+fastq[1]+" -p "+str(t)+" -o "+filehandler.get_data_path()+"outsalmon", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if res.wait() != 0:
                output, error = res.communicate()
                print(error)
                print(output)
    elif aligner == "hisat2":
        print("align with hisat2")
    
    return read_result(aligner)

def read_result(aligner):
    res = ""
    if aligner == "kallisto":
        res = pd.read_csv(filehandler.get_data_path()+"outkallisto/abundance.tsv", sep="\t")
    elif aligner == "salmon":
        res = pd.read_csv(filehandler.get_data_path()+"outsalmon/quant.sf", sep="\t")
    elif aligner == "hisat2":
        print("missing")
    return res