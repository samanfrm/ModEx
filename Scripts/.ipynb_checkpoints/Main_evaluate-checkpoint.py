#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###confusion matrix

import pandas as pd
import functions as fn
from configCHIP import output_directory
from configCHIP import input_directory
from configCHIP import result_directory
from sklearn.metrics import confusion_matrix

result_directory='Data/CHIP_split/results/'

raw_result=pd.read_csv(result_directory+"ChIPfilter-1K-950_rels_result.csv", sep="\t",header=(0),dtype=object)
raw_result.fillna('', inplace=True)
no_detected_pmid=raw_result[raw_result["find_pmid"]=='0']
detected_pmid=raw_result[raw_result["find_pmid"]!='0']
mode_detected=detected_pmid[detected_pmid["score"]!='0']
positive_mode=detected_pmid[detected_pmid["score"]>'0']
negative_mode=detected_pmid[detected_pmid["score"]<'0']

print('Number of all interactions: ',len(raw_result))
print('Number of undetected interactions: ', len(no_detected_pmid))
print('Number of detected interactions: ' ,len(detected_pmid))
print('Number of detected interactions with mode of regulation: ',
len(mode_detected))
print('Number of positive mode: ',len(positive_mode))
print('Number of negative mode: ',len(negative_mode))





