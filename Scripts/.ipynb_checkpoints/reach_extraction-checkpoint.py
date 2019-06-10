#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from Bio import Entrez
from Rules_Class import Rules
import functions as fn
from sklearn.metrics import confusion_matrix
import sys
from config import output_directory
from config import input_directory
from config import result_directory
#import config
#part=sys.argv[1]
#part='part0'
#trrust_raw=pd.read_csv(output_directory+part+".csv",sep='\t',header=(0),dtype=object)
#
#trrust_raw_result=pd.DataFrame(columns=['SrcID','SrcName','TrgID','TrgName','Mode','PMID','output','flag'])
#
#
#Positive=[]
#[Positive.append(line.strip().upper()) for line in open(input_directory+"Positive.txt")]
#Negative=[]
#[Negative.append(line.strip().upper()) for line in open(input_directory+"Negative.txt")]
#
#genes_ents=input_directory + "ALL_Human_Genes_Info.csv" #NCBI
#genes=pd.read_csv(genes_ents,sep=',',header=(0))
#genes.fillna('', inplace=True)
#
#lookup_ids=pd.read_csv(input_directory+"ncbi_id_lookup.csv",sep='\t',header=(0))
#
#for row in trrust_raw.itertuples():
#    #try:
#        query_id=[int(row.SrcID),int(row.TrgID)] #NCBI ID MUST [TF,Target]
#        query_genes=[row.SrcName,row.TrgName] #Symbol MUST [TF,Target]
#        print(query_genes)
#        PubIDs =row.PMID.split(';')
#        query_genes, query_id,single_gene,single_id =fn.make_query(genes,lookup_ids,query_genes,query_id)
#        myterm=fn.term_maker(single_gene,genes)
#        Entrez.email="saman.farahmand001@umb.edu"
#
##        try:
##            handle=Entrez.esearch(db="pubmed", term=myterm, retmax=100000000)
##            record=Entrez.read(handle)
##        except:
##            print("Enterz Fetch Error!")
#            #trrust_raw_result=trrust_raw_result.append({'SrcID':row.SrcID,'SrcName':row.SrcName,'TrgID':row.TrgID,'TrgName':row.TrgName,'Mode':row.Mode,'PMID':row.PMID,'output':'error','rank':float('nan'),'status':'Enterz Fetch Error'},ignore_index=True)
#  #          continue
#        mode=''
#        flag='0'
#        if(len(PubIDs)>0):
#            sum_ranks=[]
#            for PubID in PubIDs:
#                status=''
#                ranks=[]
#                annot_df=pd.DataFrame(columns=['type','id','text','offset','length'])
#                try:
#                    annot_df, abstract=fn.pubtator_annot(annot_df,PubID)
#                except:
#                    abstract=fn.ret_abstract(PubID)
#                    if(abstract=='?'):
#                        print("PMID=["+PubID+"] does not exist any more!")
#                        continue # remove it from the output results in TRRUST
#                    else:
#                        status+="PMID=["+PubID+"] PubTator Response is not readable, Try to annotate manually..."
#                        print(status)
#                    #add surface similarity
#                try:
#                    beCAS_lookup_full=fn.beCAS_lookup(PubID,query_id)
#                    beCAS_lookup=beCAS_lookup_full[['type','id','text','offset','length']]
#                    annot_df=pd.concat([annot_df,beCAS_lookup], ignore_index=True)
#                except:
#                    print("beCAS Server error...!")
#
#                lookup_results=fn.lookup_annot(abstract,query_genes,query_id,lookup_ids)
#                annot_df=annot_df.append(lookup_results)
#                surface_annot=fn.surface_similarity(abstract, genes, query_genes, query_id,lookup_ids,single_id)
#                annot_df=annot_df.append(surface_annot)
#                annot_df=annot_df.drop_duplicates(subset=['id','offset'])
#
#                annot_df=fn.multiple_taxonomy(annot_df, query_id)
#
#
#                annot_df=annot_df.reset_index(drop=True)
#                candidate_sentences, covered=fn.candidate_sentence(annot_df,abstract,query_id)
#                if(len(candidate_sentences.index)==0):
#                    status+="PMID=["+PubID+"] No co-existed sentences found in the abstract...!"
#                    print(status)
#                    continue
#                flag,mode=fn.REACH_apiv2(candidate_sentences,query_genes,query_id,Positive,Negative)
#                print('flag value for NER detection: ', flag)
#                print('extracted mode is: ',mode)
#                if(mode is not None):
#                    break
#
#
#            trrust_raw_result=trrust_raw_result.append({'SrcID':row.SrcID,'SrcName':row.SrcName,'TrgID':row.TrgID,'TrgName':row.TrgName,'Mode':row.Mode,'PMID':row.PMID,'output':mode,'flag':flag},ignore_index=True)
#        else:
#            print("Not found any PMIDs for this interaction...!")
#            trrust_raw_result=trrust_raw_result.append({'SrcID':row.SrcID,'SrcName':row.SrcName,'TrgID':row.TrgID,'TrgName':row.TrgName,'Mode':row.Mode,'PMID':row.PMID,'output':'error','flag':float('nan')},ignore_index=True)
#    #except:
#        #print("general Error!!")
#        #continue
#
#trrust_raw_result.to_csv(result_directory+part+"-REACH-result.csv",sep = '\t')


###confusion matrix
trrust_raw_result=pd.read_csv(result_directory+"REACH_trrust_raw_result.csv", sep="\t",header=(0),dtype=object)
trrust_raw_result.fillna('', inplace=True)
trrust_not_detected=trrust_raw_result[(trrust_raw_result["output"]!='')]
trrust_detected=trrust_raw_result[(trrust_raw_result["output"]!='')]
print('Number of all interactions:', len(trrust_raw_result))
print('Number of dettected: ',len(trrust_detected))
test=trrust_detected["Mode"].tolist()
pred=trrust_detected["output"].tolist()
cnf_matrix = confusion_matrix(test, pred)
print(cnf_matrix)
print("accuracy considering identified interactions: " + str((cnf_matrix[0][0]+cnf_matrix[1][1])/len(trrust_detected.index)*100))
print("accuracy considering all interactions: " + str((cnf_matrix[0][0]+cnf_matrix[1][1])/len(trrust_raw_result.index)*100))







