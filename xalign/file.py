import numpy as np
import gzip
import matplotlib.pyplot as plt
import pandas as pd
import requests, sys
import urllib.request
import os
import shutil

def download_file(url: str, target, overwrite: bool = False, verbose: bool = False) -> str:
    if (not os.path.exists(get_data_path()+target) or overwrite):
        if verbose:
            print("Download file")
        urllib.request.urlretrieve(url, get_data_path()+target)
    else:
        if verbose:
            print("File cached. To reload use download_file(\""+url+"\", \""+target+"\", overwrite=True) instead.")
    return get_data_path()+target

def gunzip(source: str, target, overwrite: bool = False, remove_source: bool = True, verbose: bool = False):
    if (not os.path.exists(get_data_path()+target) or overwrite):
        if verbose:
            print("gunzip file")
        with gzip.open(get_data_path()+source, 'rb') as fr:
            with open(get_data_path()+target, 'wb') as fw:
                shutil.copyfileobj(fr, fw)
        if remove_source:
            os.unlink(source)
    else:
        if verbose:
            print("File cached. To reload use gunzip(\""+source+"\", \""+target+"\", overwrite=True) instead.")
    return get_data_path()+target

def get_data_path() -> str:
    path = os.path.join(
        os.path.dirname(__file__),
        'data/'
    )
    return(path)

def concat(file_1: str, file_2: str, verbose: bool = False):
    if verbose:
            print("Concat references coding/non-coding")
    f = gzip.open(get_data_path()+file_2)
    file_content = f.read()
    f.close()
    f = gzip.open(get_data_path()+file_1, 'a', 9)
    f.write(file_content)
    f.close()
