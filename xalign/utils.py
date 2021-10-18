import gzip
import numpy as np
import os
import glob

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

def overlap(id1, id2):
    for i in id2:
        if i in id1:
            return True
    return False

def find_match(file, n=1000):
    filepath = os.path.dirname(file)
    print(filepath)
    files = glob.glob(filepath+"/*.fastq*")
    files.remove(file)
    id1 = get_readids(file, n)
    id1 = set(id1)
    for f in files:
        id2 = get_readids(f, n)
        if overlap(id1, id2):
            return f
    return ""

def file_pairs(filepath, n=1000):
    files = sorted(glob.glob(filepath+"/*.fastq*"))
    pairs = []
    for f in files:
        fm = find_match(f, n)
        if fm in files:
            files.remove(fm)
        pairs.append((f,fm))
    return pairs
        