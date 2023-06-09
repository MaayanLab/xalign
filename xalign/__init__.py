import pandas as pd
import sys
import os
import sys
import platform
import tarfile
import subprocess
from tqdm import tqdm
import re

import xalign.file as filehandler
import xalign.ensembl as ensembl
import xalign.sra as sra
from xalign.ensembl import retrieve_ensembl_organisms, organism_display_to_name
from xalign.utils import file_pairs

def build_index(aligner: str, species: str, release=None, noncoding=False, overwrite=False, verbose=False):
    organisms = ensembl.retrieve_ensembl_organisms(str(release))
    if species in organisms:
        if verbose:
            print("Download fastq.gz")
            print(filehandler.get_data_path())
        filehandler.download_file(organisms[species][2], species+"."+str(release)+".fastq.gz", overwrite=overwrite, verbose=False)
        if noncoding:
            filehandler.download_file(organisms[species][4], species+"."+str(release)+".nc.fastq.gz", overwrite=overwrite, verbose=False)
            if aligner == "kallisto":
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx")) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
            elif aligner == "salmon":
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species)) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
    else:
        print("Species not found in the Ensembl database")
        sys.exit(0)
    osys = platform.system().lower()
    download_aligner(aligner, osys)
    os.makedirs(filehandler.get_data_path()+"/index", exist_ok=True)
    os.makedirs(filehandler.get_data_path()+"/index/"+str(release), exist_ok=True)
    if aligner == "kallisto":
        if (not os.path.exists(filehandler.get_data_path()+"index/"+str(str(release))+"/kallisto_"+species+".idx")) or overwrite:
            if verbose:
                print("Build kallisto index for "+species)
            os.system(filehandler.get_data_path()+"kallisto/kallisto index -i "+filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx "+filehandler.get_data_path()+species+"."+str(release)+".fastq.gz")
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")
    elif aligner == "salmon":
        if (not os.path.exists(filehandler.get_data_path()+"index/salmon_"+species)) or overwrite:
            if verbose:
                print("Build salmon index for "+species)
                print(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin index -i "+filehandler.get_data_path()+"index/"+str(str(release))+"/salmon_"+species+".idx -t "+filehandler.get_data_path()+species+"."+str(release)+".fastq.gz")
            os.system(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon index -i "+filehandler.get_data_path()+"index/"+str(str(release))+"/salmon_"+species+" -t "+filehandler.get_data_path()+species+"."+str(release)+".fastq.gz")
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")

def download_aligner(aligner, osys, verbose=False):
    if verbose:
        print(osys)

    if aligner == "salmon":
        if osys == "linux":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extractall(filehandler.get_data_path())
            file.close()
        if osys == "darwin":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extractall(filehandler.get_data_path())
            file.close()
        else:
            if verbose:
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

def align_fastq(species, fastq, aligner="kallisto", t=1, release=None, noncoding=False, overwrite=False, verbose=False):
    if isinstance(fastq, str):
        fastq = [fastq]
    if not release:
        release = list(ensembl.retrieve_ensembl_organisms(release).items())[0][1][5]
    build_index(aligner, species, release=release, noncoding=noncoding, overwrite=overwrite, verbose=verbose)
    if aligner == "kallisto":
        if len(fastq) == 1:
            if verbose:
                print("Align with kallisto (single strand).")
            res = subprocess.Popen(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto --single -l 200 -s 20 "+fastq[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        else:
            if verbose:
                print("Align with kallisto (paired).")
            res = subprocess.Popen(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto "+fastq[0]+" "+fastq[1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if res.wait() != 0:
            output, error = res.communicate()
            if verbose:
                print(output)
            print(error)
    elif aligner == "salmon":
        if len(fastq) == 1:
            print("Align with salmon (single).")
            res = subprocess.Popen(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species+" -l A -r "+fastq[0]+" -p "+str(t)+" --validateMappings -o "+filehandler.get_data_path()+"outsalmon", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if res.wait() != 0:
                output, error = res.communicate()
                print(error)
                if verbose:
                    print(output)
        else:
            print("Align with salmon (paired).")
            res = subprocess.Popen(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species+" -l A -1 "+fastq[0]+" -2 "+fastq[1]+" -p "+str(t)+" -o "+filehandler.get_data_path()+"outsalmon", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if res.wait() != 0:
                output, error = res.communicate()
                print(error)
                if verbose:
                    print(output)
    elif aligner == "hisat2":
        print("align with hisat2")
    
    return read_result(aligner)

def align_folder(species, folder, aligner="kallisto", t=1, release=None, identifier="symbol", noncoding=False, overwrite=False, verbose=False):
    fastq_files = file_pairs(folder)
    gene_counts = []
    transcript_counts = []
    pbar = tqdm(total=len(fastq_files))
    pbar.set_description("Aligning samples")
    sample_names = []
    for fq in fastq_files:
        if fq[0] == "" or fq[1] == "":
            fq.remove("")
            sample_names.append(fq[0])
        else:
            bnames = [os.path.basename(x) for x in fq]
            sample_names.append(re.sub(r'_$','',os.path.commonprefix(bnames)))
        res = align_fastq(species, fq, aligner=aligner, t=t, release=release, noncoding=noncoding, overwrite=overwrite, verbose=verbose)
        transcript_counts.append(list(res.loc[:,"reads"].round()))
        res_gene = ensembl.agg_gene_counts(res, species, identifier=identifier, overwrite=overwrite)
        gene_counts.append(list(res_gene.loc[:,"counts"].round()))
        overwrite=False
        pbar.update(1)
    transcript_counts = pd.DataFrame(transcript_counts, columns=res.iloc[:,0], index=sample_names).T.astype("int")
    gene_counts = pd.DataFrame(gene_counts, columns=res_gene.iloc[:,0], index=sample_names).T.astype("int")
    return (gene_counts, transcript_counts)

def read_result(aligner):
    res = ""
    if aligner == "kallisto":
        res = pd.read_csv(filehandler.get_data_path()+"outkallisto/abundance.tsv", sep="\t")
        res = res.loc[:,["target_id", "est_counts", "tpm"]]
        res.columns = ["transcript", "reads", "tpm"]
    elif aligner == "salmon":
        res = pd.read_csv(filehandler.get_data_path()+"outsalmon/quant.sf", sep="\t")
        res = res.loc[:,["Name", "NumReads", "TPM"]]
        res.columns = ["transcript", "reads", "tpm"]
    elif aligner == "hisat2":
        print("missing")
    return res