import subprocess
import multiprocessing
from tqdm import tqdm
import tarfile
import os

import xalign.file as filehandler
import xalign.utils

def load_sra_star(args):
    load_sra(*args)

def load_sra(sample, output):

    if not os.path.exists(filehandler.get_data_path()+"fasterq-dump-mac"):
        file = tarfile.open(filehandler.get_data_path()+"fasterq.tar.gz")
        file.extractall(filehandler.get_data_path())
        file.close()

    if xalign.utils.current_os() == "linux":
        res = subprocess.Popen(filehandler.get_data_path()+"fasterq/fasterq-dump-ubuntu -f --mem 2G --split-3 --threads 2 --skip-technical -O "+output+" "+sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    elif xalign.utils.current_os() == "mac":
        res = subprocess.Popen(filehandler.get_data_path()+"fasterq/fasterq-dump-mac -f --mem 2G --split-3 --threads 2 --skip-technical -O "+output+" "+sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if res.wait() != 0:
        output, error = res.communicate()
        print(error)
        print(output)
    return 1

def load_sras(samples, output, t=4):
    with multiprocessing.Pool(t) as pool:
        args = [(sample, output) for sample in samples]
        res = list(tqdm(pool.imap(load_sra_star, args), desc="Downloading", total=len(args)))
    return res
