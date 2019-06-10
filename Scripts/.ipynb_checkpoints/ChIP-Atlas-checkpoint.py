#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os


directory="~/Dropbox/UMB/Dr-Koroush/Code/Data/"
chip_ents=directory + "ChIP.csv"
chip_rels=directory+"ChIP-rels.csv"

rels=pd.read_csv(chip_rels,sep='\t',header=(0))
ents=pd.read_csv(chip_ents,sep='\t',header=(0))

a=ents[ents['uid'].isin([2,766])]['name']
print(a)
