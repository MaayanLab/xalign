import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import requests, sys
import urllib.request
import os

def download_file(url: str, target, overwrite: bool = False, verbose: bool = False) -> str:
    if (not os.path.exists(get_data_path()+target) or overwrite):
        if verbose:
            print("Download file")
        urllib.request.urlretrieve(url, get_data_path()+target)
    else:
        if verbose:
            print("File cached. To reload use download_file(\""+url+"\", \""+target+"\", overwrite=True) instead.")
    return get_data_path()+target

def get_data_path() -> str:
    path = os.path.join(
        os.path.dirname(__file__),
        'data/'
    )
    return(path)
