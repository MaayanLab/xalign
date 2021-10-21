# xAlign: Hassle-free transcript quantification

xAlign is an efficient python package to align FASTQ files against any Ensembl reference genomes. The currently supported alignment algorithms are `kallisto` (https://pachterlab.github.io/kallisto/) and `Salmon` (https://salmon.readthedocs.io/en/latest/salmon.html). The package contains modules for Ensemble ID mapping to gene symbols via the `mygene.info` python package and SRA download capabilities. When using this package please cite the corresponding alignment algorithm.

## Installation

```
pip3 install git+https://github.com/MaayanLab/xalign.git
```

## Requirements

The alignment algorithms require a minimum of around 5GB of memory to run. When downloading SRA files, make sure that there is sufficient available disk space. `xalign` is currently only working on `Linux` operating systems.

## Usage

The recommended usage is `xalign.align_folder()` if there are multiple FASTQ files. These FASTQ files can be aligned one by one, and gene level counts can be aggregated using the function `xalign.ensembl.agg_gene_counts()`

### Align a single FASTQ file in single-read mode

To align a single RNA-seq file we first download an example SRA file and save it in the folder `data/example_1` relative to the working directory. The function `xalign.align_fastq()` will generate the required cDNA index from the Ensembl reference genome when the index is not already built. `result` is a dataframe with transcript IDs, gene counts, and TPM.

When the alignment is run against a new species, the initial setup will take a few minutes to complete because building a new index and creating gene mapping files are required.

```python

import xalign

xalign.sra.load_sras(["SRR14457464"], "data/example_1")

result = xalign.align_fastq("homo_sapiens", "data/example_1/SRR14457464.fastq", t=8)

```

### Align a single FASTQ file in paired-end mode

To align a single RNA-seq file in paired-end mode we first download an example SRA file and save it in folder `data/example_2` relative to the working directory. If the SRA file is a paired-end sample, two files will be generated with the two suffixes `_1` and `_2`. The function `xalign.align_fastq()` will generate the required cDNA index from the Ensembl reference genome when the index is not already built. `result` is a dataframe with transcript IDs, gene counts, and TPM.

When the alignment is run against a new species, the initial setup will take a couple of minutes to built the index and to create the gene mapping files.

```python

import xalign

# the sample is paired-end and will result in two files (SRR15972519_1.fastq, SRR15972519_2.fastq)
xalign.sra.load_sras(["SRR15972519"], "data/example_2")

result = xalign.align_fastq("homo_sapiens", ["data/example_2/SRR15972519_1.fastq", "data/example_2/SRR15972519_2.fastq"], t=8)

```

### Align FASTQ files in a directory

`xalign` can automatically align all files in a given folder, instead of calling `xalign.align_fastq()` multiple times. In this case `xalign.align_folder()` will automatically detect whether the folder contains paired- or single-end samples and group the samples accordingly without manual input. The output will be two dataframes. `gene_count` will contain gene level counts that can be aggregated for different gene identifiers (symbol:default, ensembl_id, entrezgene_id). Transcripts that can not be mapped to corresponding identifiers are discarded. `transcript_count` contains the read counts at transcript level.

```python

import xalign

# this will download multiple GB of samples
xalign.sra.load_sras(["SRR15972519", "SRR15972520", "SRR15972521"], "data/example_3")

gene_count, transcript_count = xalign.align_folder("homo_sapiens", "data/example_3", t=8, overwrite=False)

```

### Mapping transcript counts to gene-level counts

When FASTQ files are aligned individually using `xalign.align_fastq()` the output is in transcript-level. To aggregate counts to gene-level the function `xalign.ensembl.agg_gene_counts()` can be used.

```python

import xalign

xalign.sra.load_sras(["SRR14457464"], "data/example_4")

result = xalign.align_fastq("homo_sapiens", "data/example_4/SRR15972519.fastq", t=8)

# identifier can be symbol/ensembl_id/entrezgene_id
gene_counts = xalign.ensembl.agg_gene_counts(result, "homo_sapiens", identifier="symbol")

```
