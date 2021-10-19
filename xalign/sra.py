import subprocess
import multiprocessing
from tqdm import tqdm
import gzip
import os

import xalign.file as filehandler

def gunzip(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)

def load_sra_star(args):
    load_sra(*args)

def load_sra(sample, output):

    if not os.path.exists(filehandler.get_data_path()+"fasterq-dump"):
        gunzip(filehandler.get_data_path()+"fasterq-dump.gz", filehandler.get_data_path()+"fasterq-dump")

    res = subprocess.Popen(filehandler.get_data_path()+"fasterq-dump -f --mem 2G --split-3 --threads 2 --skip-technical -O "+output+" "+sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if res.wait() != 0:
        output, error = res.communicate()
        print(error)
        print(output)
    return 1

def load_sras(samples, output, t=4):
    with multiprocessing.Pool(t) as pool:
        args = [(sample, output) for sample in samples]
        print(args)
        res = list(tqdm(pool.imap(load_sra_star, args), desc="Downloading", total=len(args)))
    return res
