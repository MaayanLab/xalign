import subprocess

import xalign.file as filehandler

def load_sra(sample, output):
    res = subprocess.Popen(filehandler.get_data_path()+"fasterq-dump --split-files -O "+output+" "+sample, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if res.wait() != 0:
        output, error = res.communicate()
        print(error)
        print(output)