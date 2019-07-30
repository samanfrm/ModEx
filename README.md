# ModEx
A text mining system for extracting mode of regulation of Transcription Factor-gene regulatory interaction

## Overview

Deciphering the network of TF-target interactions with information on mode of regulation (activation vs. repression) is an important step toward understanding the regulatory pathways that underlie complex traits. There are many experimental, computational, and manually curated databases of TF-gene interactions. In particular, high-throughput ChIP-seq datasets provide a large-scale map or transcriptional regulatory interactions.

However, these interactions are not annotated with information on context and mode of regulation. Such information is crucial to gain a global picture of gene regulatory mechanisms and can aid in developing machine learning models for applications such as biomarker discovery, prediction of response to therapy, and precision medicine. 

we introduce a text-mining system, to annotate ChIP-seq derived interaction with such meta data through mining PubMed articles.

## Dependencies and Installation
ModEx can be installed using the GitHub repository. All of the dependencies will be installed via setup.py script.
```
git clone https://github.com/samanfrm/ModEx
cd gbnet
python3 setup.py install --user
cd ..
```
## Import libraries
We need to import required libraries into the script:
```python
import pandas as pd
import functions as fn
import os
```
## Run the system for a sample query interaction
First, the paths to necessary files and dictionaries must be defined based on relative path from script directory:

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
```
Then, we need to create the query variables and assign them with the transcription factor and target genes entrez IDs respectively:

```python
# [TF_ID, Target_ID]
query_id=[26574,4609]
```

Next, we need to set the binding port to [Stanford CoreNLP](https://stanfordnlp.github.io/CoreNLP/):

```python
parser_port="8000"
```
Also, optional values for the MeSH term and email address should be defined:

```python
mesh='humans'
email='example@mail.com'
```
Finally, we can run the test mining system to annotate the query interaction as well as associated evidence and citations:

```python
res=fn.modex(query_id,parser_port,Positive,Negative,lookup_ids,genes,mesh,email)
```
The result is a dataframe including mode of regulation and all of the associated citations and evidence sentences for the annotation:

| src_entrez  |  trg_entrez | srcname  | trgname  |  mode     | score  | evi_pmid        | evi_sent                  |
|-------------|-------------|----------|----------|-----------|--------|-----------------|---------------------------|
|  26574      |  4609       | AATF     |  MYC     |  positive | 4      | 20924650;2054...| [20924650]WE HAVE UNAMB...|


## Citation

Saman Farahmand, Todd Riley, Kourosh Zarringhalam, "ModEx: A text mining system for extracting mode of regulation of Transcription Factor-gene regulatory interaction", BioRxiv, 2019
