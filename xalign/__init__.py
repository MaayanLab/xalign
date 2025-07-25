import pandas as pd
import sys
import os
import sys
import platform
import tarfile
import zipfile
import subprocess
import typing as t
from tqdm import tqdm
import re
import shutil

import xalign.file as filehandler
import xalign.ensembl as ensembl
import xalign.sra as sra
from xalign.ensembl import retrieve_ensembl_organisms, organism_display_to_name
from xalign.utils import file_pairs

Aligner = t.Literal['kallisto', 'salmon', 'hisat2', 'star']

def build_index(aligner: Aligner, species: str, release=None, noncoding=False, overwrite=False, verbose=False, t=1):
    organisms = ensembl.retrieve_ensembl_organisms(str(release))
    if species in organisms:
        if verbose:
            print("Download fastq.gz")
            print(filehandler.get_data_path())

        filehandler.download_file(organisms[species]['cdna_url'], species+"."+str(release)+".fastq.gz", overwrite=overwrite, verbose=verbose)
        if noncoding:
            filehandler.download_file(organisms[species]['ncdna_url'], species+"."+str(release)+".nc.fastq.gz", overwrite=overwrite, verbose=verbose)
            if aligner == "kallisto":
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx")) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
            elif aligner == "salmon":
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species)) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
            elif aligner == "hisat2":
                # TODO: do we need to do anything else for hisat2?
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/hisat2_"+species)) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
            elif aligner == "star":
                if (not os.path.exists(filehandler.get_data_path()+"index/"+str(release)+"/star_"+species)) or overwrite:
                    filehandler.concat(species+"."+str(release)+".fastq.gz", species+"."+str(release)+".nc.fastq.gz", verbose=verbose)
            else:
                raise NotImplementedError(aligner)

        if aligner in "star":
            filehandler.download_file(organisms[species]['gtf_url'], species+"."+str(release)+".gtf.gz", overwrite=overwrite, verbose=verbose)
            filehandler.download_file(organisms[species]['primary_assembly_fa_url'], species+"."+str(release)+".fa.gz", overwrite=overwrite, verbose=verbose)
            filehandler.gunzip(species+"."+str(release)+".fa.gz", species+"."+str(release)+".fa", overwrite=overwrite, verbose=verbose)
            filehandler.gunzip(species+"."+str(release)+".gtf.gz", species+"."+str(release)+".gtf", overwrite=overwrite, verbose=verbose)
    else:
        print("Species not found in the Ensembl database")
        sys.exit(0)
    osys = platform.system().lower()
    download_aligner(aligner, osys)
    os.makedirs(filehandler.get_data_path()+"/index", exist_ok=True)
    os.makedirs(filehandler.get_data_path()+"/index/"+str(release), exist_ok=True)
    if aligner == "kallisto":
        if (not os.path.exists(filehandler.get_data_path()+"index/"+str(str(release))+"/kallisto_"+species+".idx")) or overwrite:
            args = [
                filehandler.get_data_path()+"kallisto/kallisto", "index",
                "-i", filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx",
                filehandler.get_data_path()+species+"."+str(release)+".fastq.gz",
            ]
            if verbose:
                print("Build kallisto index for "+species)
                print(*args)
            subprocess.run(args, stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")
    elif aligner == "salmon":
        if (not os.path.exists(filehandler.get_data_path()+"index/salmon_"+species)) or overwrite:
            args = [
                filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon", "index",
                "-p", str(t),
                "-i", filehandler.get_data_path()+"index/"+str(str(release))+"/salmon_"+species,
                "-t", filehandler.get_data_path()+species+"."+str(release)+".fastq.gz",
            ]
            if verbose:
                print("Build salmon index for "+species)
                print(*args)
            subprocess.run(args, stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")
    elif aligner == "hisat2":
        if (not os.path.exists(filehandler.get_data_path()+"index/hisat2_"+species)) or overwrite:
            args = [
                shutil.which("hisat2-build"),
                "-p", str(t),
                # TODO: where do we get GENOME/FA
                filehandler.get_data_path()+species+"."+str(release)+".fa",
                filehandler.get_data_path()+"index/"+str(str(release))+"/hisat2_"+species,
            ]
            if verbose:
                print("Build hisat2 index for "+species)
                print(*args)
            subprocess.run(args, stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")
    elif aligner == "star":
        if (not os.path.exists(filehandler.get_data_path()+"index/star_"+species)) or overwrite:
            args = [
                filehandler.get_data_path()+"STAR",
                "--runMode", "genomeGenerate",
                "--genomeDir", filehandler.get_data_path()+"index/"+str(str(release))+"/star_"+species,
                "--genomeFastaFiles", filehandler.get_data_path()+species+"."+str(release)+".fa",
                "--sjdbGTFfile", filehandler.get_data_path()+species+"."+str(release)+".gtf",
                "--runThreadN", str(t),
                "--genomeSAsparseD", "2",
                "--genomeSAindexNbases", "13",
            ]
            if verbose:
                print("Build star index for "+species)
                print(*args)
            subprocess.run(args, stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            if verbose:
                print("Index already exists. Use overwrite to rebuild.")
    else:
        raise NotImplementedError(aligner)

def download_aligner(aligner: Aligner, osys, verbose=False):
    if verbose:
        print(osys)

    if aligner == "salmon":
        if osys == "linux":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extractall(filehandler.get_data_path())
            file.close()
        elif osys == "darwin":
            url = "https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz"
            filepath = filehandler.download_file(url, "salmon.tar.gz")
            file = tarfile.open(filepath)
            file.extractall(filehandler.get_data_path())
            file.close()
        else:
            if verbose:
                print("Salmon not supported by this package for this operating system.")

    elif aligner == "kallisto":
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
        assert shutil.which('hisat2'), 'Please install hisat yourself'
        assert shutil.which("hisat2-build"), 'Please install hisat-build yourself'
        assert shutil.which('featureCounts'), 'Please install featureCounts yourself'

    elif aligner == "star":
        if shutil.which('STAR'):
            os.link(shutil.which('STAR'), filehandler.get_data_path()+"STAR")
        elif os.path.exists(filehandler.get_data_path()+"STAR"):
            pass
        elif osys == "windows":
            assert shutil.which('STAR'), 'Please install star yourself'
        elif osys == "linux":
            url = "https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip"
            filepath = filehandler.download_file(url, "STAR.zip")
            file = zipfile.ZipFile(filepath)
            file.extract('STAR_2.7.11b/Linux_x86_64_static/STAR', filehandler.get_data_path())
            file.close()
            os.link(filehandler.get_data_path()+'STAR_2.7.11b/Linux_x86_64_static/STAR', filehandler.get_data_path()+"STAR")
            os.chmod(filehandler.get_data_path()+"STAR", 0o700)
        else: #mac
            url = "https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip"
            filepath = filehandler.download_file(url, "STAR.zip")
            file = zipfile.ZipFile(filepath)
            file.extract('STAR_2.7.11b/MacOSX_x86_64/STAR', filehandler.get_data_path())
            file.close()
            os.link(filehandler.get_data_path()+'STAR_2.7.11b/MacOSX_x86_64/STAR', filehandler.get_data_path()+"STAR")
            os.chmod(filehandler.get_data_path()+"STAR", 0o700)

    else:
        raise NotImplementedError(aligner)

def align_fastq(species, fastq, aligner: Aligner="kallisto", t=1, release=None, noncoding=False, overwrite=False, verbose=False):
    if isinstance(fastq, str):
        fastq = [fastq]
    if not release:
        release = list(ensembl.retrieve_ensembl_organisms(release).items())[0][1]['release']
    build_index(aligner, species, release=release, noncoding=noncoding, overwrite=overwrite, verbose=verbose)
    if aligner == "kallisto":
        if len(fastq) == 1:
            if verbose:
                print("Align with kallisto (single strand).")
            subprocess.run(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto --single -l 200 -s 20 "+fastq[0], stdout=sys.stdout if verbose else None, stderr=sys.stderr, shell=True)
        else:
            if verbose:
                print("Align with kallisto (paired).")
            subprocess.run(filehandler.get_data_path()+"kallisto/kallisto quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/kallisto_"+species+".idx -t "+str(t)+" -o "+filehandler.get_data_path()+"outkallisto "+fastq[0]+" "+fastq[1], stdout=sys.stdout if verbose else None, stderr=sys.stderr, shell=True)
    elif aligner == "salmon":
        if len(fastq) == 1:
            print("Align with salmon (single).")
            subprocess.run(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species+" -l A -r "+fastq[0]+" -p "+str(t)+" --validateMappings -o "+filehandler.get_data_path()+"outsalmon", stdout=sys.stdout if verbose else None, stderr=sys.stderr, shell=True, check=True)
        else:
            print("Align with salmon (paired).")
            subprocess.run(filehandler.get_data_path()+"salmon-1.5.2_linux_x86_64/bin/salmon quant -i "+filehandler.get_data_path()+"index/"+str(release)+"/salmon_"+species+" -l A -1 "+fastq[0]+" -2 "+fastq[1]+" -p "+str(t)+" -o "+filehandler.get_data_path()+"outsalmon", stdout=sys.stdout if verbose else None, stderr=sys.stderr, shell=True, check=True)
    elif aligner == "hisat2":
        if len(fastq) == 1:
            print("Align with hisat2 (single).")
            subprocess.run([
                shutil.which("hisat2"),
                "-x", filehandler.get_data_path()+"index/"+str(release)+"/hisat2_"+species,
                "-U", fastq[0],
                "-p", str(t),
                "-S", filehandler.get_data_path()+"outhisat2/out.sam",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
            subprocess.run([
                shutil.which("featureCounts"),
                "-T", str(t),
                # TODO: what's up with this file
                "-a", filehandler.get_data_path()+"index/"+str(release)+"/hisat2_"+species+".gtf",
                "-o", filehandler.get_data_path()+"outhisat2/out.tsv",
                "-S", filehandler.get_data_path()+"outhisat2/out.sam",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            print("Align with hisat2 (paired).")
            subprocess.run([
                shutil.which("hisat2"),
                "-x", filehandler.get_data_path()+"index/"+str(release)+"/hisat2_"+species,
                "-1", fastq[0],
                "-2", fastq[1],
                "-p", str(t),
                "-S", filehandler.get_data_path()+"outhisat2/out.sam",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
            subprocess.run([
                shutil.which("featureCounts"),
                "-T", str(t),
                # TODO: what's up with this file
                "-a", filehandler.get_data_path()+"index/"+str(release)+"/hisat2_"+species+".gtf",
                "-o", filehandler.get_data_path()+"outhisat2/out.tsv",
                "-S", filehandler.get_data_path()+"outhisat2/out.sam",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
    elif aligner == "star":
        if len(fastq) == 1:
            print("Align with star (single).")
            subprocess.run([
                filehandler.get_data_path()+"STAR",
                "--genomeDir", filehandler.get_data_path()+"index/"+str(release)+"/star_"+species,
                "--limitBAMsortRAM", "10000000000",
                "--runThreadN", str(t),
                "--outSAMstrandField", "intronMotif",
                "--outFilterIntronMotifs", "RemoveNoncanonical",
                "--outFileNamePrefix", filehandler.get_data_path()+"outstar",
                "--readFilesIn", fastq[0],
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outReadsUnmapped", "Fastx",
                "--outSAMmode", "Full",
                "--quantMode", "GeneCounts",
                "--limitIObufferSize", "50000000", "50000000",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
        else:
            print("Align with star (paired).")
            subprocess.run([
                filehandler.get_data_path()+"STAR",
                "--genomeDir", filehandler.get_data_path()+"index/"+str(release)+"/star_"+species,
                "--limitBAMsortRAM", "10000000000",
                "--runThreadN", str(t),
                "--outSAMstrandField", "intronMotif",
                "--outFilterIntronMotifs", "RemoveNoncanonical",
                "--outFileNamePrefix", filehandler.get_data_path()+"outstar",
                "--readFilesIn", fastq[0], fastq[1],
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outReadsUnmapped", "Fastx",
                "--outSAMmode", "Full",
                "--quantMode", "GeneCounts",
                "--limitIObufferSize", "50000000", "50000000",
            ], stdout=sys.stdout if verbose else None, stderr=sys.stderr, check=True)
    else:
        raise NotImplementedError(aligner)
    
    return read_result(aligner)

def align_folder(species, folder, aligner: Aligner="kallisto", t=1, release=None, identifier="symbol", noncoding=False, overwrite=False, verbose=False):
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

def read_result(aligner: Aligner):
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
        res = pd.read_csv(filehandler.get_data_path()+"outhisat2/out.tsv", sep="\t")
        # TODO
        # res = res.loc[:,["Name", "NumReads", "TPM"]]
        # res.columns = ["transcript", "reads", "tpm"]
    elif aligner == "star":
        res = pd.read_csv(filehandler.get_data_path()+"outstarReadsPerGene.out.tab", sep="\t", skiprows=4, header=None)
        res.columns = ["transcript", "reads", "stranded_reads_1", "stranded_reads_2"]
        # TODO: get tpm for consistency
        res = res.loc[:,["transcript", "reads"]]
    else:
        raise NotImplementedError(aligner)
    return res
