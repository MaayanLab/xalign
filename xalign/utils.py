import gzip
import numpy as np
import glob
import platform
import functools

@functools.lru_cache()
def get_readids(file, n):
    n = n*4
    head = []
    if file.endswith(".gz"):
        with gzip.open(file, 'r') as f:
            head = [next(f).decode("UTF-8") for x in range(n)]
    else:
        with open(file) as f:
            head = [next(f) for x in range(n)]
            
    lines = np.array(head)[np.arange(0, len(head), 4)]
    
    ids = []
    for line in lines:
        sp = line.split(" ")
        ids.append(sp[0].replace("@", ""))
    return sorted(ids)

def current_os():
    osys = platform.system().lower()
    if osys == "darwin":
        return "mac"
    else:
        return osys
    
def find_match(files, file, n=1000):
    files = files - {file}
    id1 = set(get_readids(file, n))
    for f in files:
        id2 = set(get_readids(f, n))
        if id1 & id2:
            return f
    return ""

def file_pairs(filepath, n=1000):
    files = set(glob.glob(filepath+"/*.fastq*") if type(filepath) == str else filepath)
    done = set()
    pairs = []
    for f in files:
        if f in done: continue
        fm = find_match(files, f, n)
        done.add(fm)
        pairs.append([f,fm])
    return pairs
