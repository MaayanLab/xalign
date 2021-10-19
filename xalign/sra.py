import subprocess
import multiprocessing
from tqdm import tqdm

import xalign.file as filehandler

def load_sra_star(args):
    load_sra(*args)

def load_sra(sample, output):
    res = subprocess.Popen(filehandler.get_data_path()+"fasterq-dump  --mem 2G --split-3 --threads 2 --skip-technical -O "+output+" "+sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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
