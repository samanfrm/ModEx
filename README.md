# ModEx
A text mining system for extracting mode of regulation of Transcription Factor-gene regulatory interaction

## Overview

Deciphering the network of TF-target interactions with information on mode of regulation (activation vs. repression) is an important step toward understanding the regulatory pathways that underlie complex traits. There are many experimental, computational, and manually curated databases of TF-gene interactions. In particular, high-throughput ChIP-seq datasets provide a large-scale map or transcriptional regulatory interactions.

However, these interactions are not annotated with information on context and mode of regulation. Such information is crucial to gain a global picture of gene regulatory mechanisms and can aid in developing machine learning models for applications such as biomarker discovery, prediction of response to therapy, and precision medicine. 

we introduce a text-mining system, to annotate ChIP-seq derived interaction with such meta data trough mining PubMed articles.

## Dependencies and Installation
ModEx can be installed using the GitHub repository. All of the dependencies will be installed via setup.py script.
```
git clone https://github.com/samanfrm/ModEx
cd gbnet
python3 setup.py install --user
cd ..
```
## Import libraries
```python
import pandas as pd
import functions as fn
import os
```
## Run the system for a query interaction
First, the paths to necessary files and dictionaries must be defined basedn on relative path from script directory:

```python
input_directory=os.path.realpath('../Data')

Positive=[]
[Positive.append(line.strip().upper()) for line in open(input_directory+"/Positive.txt")]
Negative=[]
[Negative.append(line.strip().upper()) for line in open(input_directory+"/Negative.txt")]

genes_ents=input_directory + "/ALL_Human_Genes_Info.csv"
genes=pd.read_csv(genes_ents,sep=',',header=(0))
genes.fillna('', inplace=True)

lookup_ids=pd.read_csv(input_directory+"/ncbi_id_lookup.csv",sep='\t',header=(0))

CHIP_result=pd.DataFrame(columns=['src_entrez','trg_entrez','srcname','trgname','mode','score','evi_pmid','evi_sent'])
```
