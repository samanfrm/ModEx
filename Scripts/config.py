#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd
input_directory="../Data/"
output_directory="../Data/split/"
result_directory="../Data/split/results-with_evidence/"
input_name="trrust_with_mode"
output_name="trrust_raw_result.csv"
from sklearn.metrics import confusion_matrix
if not os.path.exists(output_directory):
    os.mkdir(output_directory)
if not os.path.exists(result_directory):
    os.mkdir(result_directory)
task=''
#task='mixscore'
#task=sys.argv[1]    #'merge' # or 'merge'
#task='merge'
if(task=='split'):
    genes_ents=input_directory + input_name #NCBI
    genes=pd.read_csv(genes_ents,sep='\t',header=(0))
    num=int(sys.argv[2])

    total=len(genes.index)
    length=total//num
    remain=total%num
    start=0

    for i in range(0,num):
        temp_df=genes[start:start+length]
        temp_df.to_csv(output_directory+'part'+str(i)+".csv",sep='\t')
        start+=length

    if(remain!=0):
        temp_df=genes[start:]
        temp_df.to_csv(output_directory+'part'+str(i+1)+".csv", sep='\t')

if(task=='merge'):
    dirlist = os.listdir(result_directory) # dir is your directory path
    number_files = len(dirlist)
    result_df=pd.read_csv(result_directory+"part0-REACH-result.csv", sep="\t",header=(0),dtype=object)
    result_df.fillna('', inplace=True)
    for i in range(1,number_files):
        temp_df=pd.read_csv(result_directory+"part"+str(i)+"-REACH-result.csv", sep="\t",header=(0),dtype=object)
        temp_df.fillna('', inplace=True)
        result_df=pd.concat([result_df,temp_df],ignore_index=True)
    result_df.to_csv(result_directory+output_name,sep = '\t')

if(task=='mixscore'):
    coef_kn=2
    coef_dpt=1
    magnit=0 # yes or o for not
#    coef_kn=int(sys.argv[2])
#    coef_dpt=int(sys.argc[3])
    trrust_raw_result_kn=pd.read_csv('trrust_raw_result-nearest.csv',sep='\t',header=(0))
    trrust_raw_result_dpt=pd.read_csv('trrust_raw_result-dpt.csv',sep='\t',header=(0))
    trrust_raw_result_kn.fillna('', inplace=True)
    trrust_raw_result_dpt.fillna('', inplace=True)
    trrust_raw_result=pd.DataFrame(columns=['SrcID','SrcName','TrgID','TrgName','Mode','PMID','output','rank','status'])
    for i in range(0,len(trrust_raw_result_dpt.index)):
        dpt_rnk=trrust_raw_result_dpt.loc[i,'rank']
        kn_rank=trrust_raw_result_kn.loc[i,'rank']
        mixs=coef_dpt*dpt_rnk + coef_kn*kn_rank
        mode=''
        if(magnit==0): #check for the bigger the better
            if(mixs>0):
                mode='positive'
            if(mixs<0):
                mode='negative'
        else:
            if(abs(dpt_rnk)>abs(kn_rank) and dpt_rnk!=0 and kn_rank!=0):
                mode=trrust_raw_result_dpt.loc[i,'output']
                mixs=dpt_rnk
            elif(abs(kn_rank)>abs(dpt_rnk) and dpt_rnk!=0 and kn_rank!=0):
                mode=trrust_raw_result_kn.loc[i,'output']
                mixs=kn_rank
            else: #if one is 0 consider another one
                if(mixs>0):
                    mode='positive'
                if(mixs<0):
                    mode='negative'
        trrust_raw_result=trrust_raw_result.append({'SrcID':trrust_raw_result_dpt.loc[i,'SrcID'],'SrcName':trrust_raw_result_dpt.loc[i,'SrcName'],'TrgID':trrust_raw_result_dpt.loc[i,'TrgID'],'TrgName':trrust_raw_result_dpt.loc[i,'TrgName'],'Mode':trrust_raw_result_dpt.loc[i,'Mode'],'PMID':trrust_raw_result_dpt.loc[i,'PMID'],'output':mode,'rank':mixs,'status':str(trrust_raw_result_dpt.loc[i,'status']) +str(trrust_raw_result_kn.loc[i,'status'])},ignore_index=True)

    trrust_raw_result.to_csv('trrust_raw_result_mixed.csv',sep = '\t')

    #trrust_raw_result=pd.read_csv("trrust_raw_result_mixed.csv", sep="\t",header=(0),dtype=object)
    #trrust_raw_result.fillna('', inplace=True)
    identified_with_mode=trrust_raw_result[trrust_raw_result['output']!='']
    print('Number of detected enteries with mode of regulation: ' + str(len(identified_with_mode.index)))
    test=identified_with_mode["Mode"].tolist()
    pred=identified_with_mode["output"].tolist()
    cnf_matrix = confusion_matrix(test, pred)
    print(cnf_matrix)
    print("accuracy considering identified interactions: " + str((cnf_matrix[0][0]+cnf_matrix[1][1])/len(identified_with_mode.index)*100))
    print("accuracy considering all interactions: " + str((cnf_matrix[0][0]+cnf_matrix[1][1])/len(trrust_raw_result.index)*100))






