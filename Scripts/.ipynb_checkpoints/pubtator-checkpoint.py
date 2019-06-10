#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
from Bio import Entrez
from Rules_Class import Rules
import functions as fn
import networkx as nx


directory="/home/saman/Dropbox/UMB/Dr-Koroush/Code/Data/"
Positive=[]
[Positive.append(line.strip().upper()) for line in open(directory+"Positive.txt")]
Negative=[]
[Negative.append(line.strip().upper()) for line in open(directory+"Negative.txt")]

genes_ents=directory + "ALL_Human_Genes_Info.csv" #NCBI
genes=pd.read_csv(genes_ents,sep=',',header=(0))
genes.fillna('', inplace=True)

lookup_ids=pd.read_csv(directory+"ncbi_id_lookup.csv",sep='\t',header=(0))


query_genes=['AATF','BAX'] #Symbol MUST [TF,Target]
query_id=[26574,581]#NCBI ID MUST [TF,Target]
PubIDs =['22909821']
query_genes, query_id,single_gene,single_id =fn.make_query(genes,lookup_ids,query_genes,query_id)
myterm=fn.term_maker(single_gene,genes)


Entrez.email="saman.farahmand001@umb.edu"
handle=Entrez.esearch(db="pubmed", term=myterm, retmax=100000000)
record=Entrez.read(handle)
IDlists = record["IdList"]
if(len(PubIDs)>0):
    sum_ranks=[]
    for PubID in PubIDs:
        ranks=[]
        annot_df=pd.DataFrame(columns=['type','id','text','offset','length'])
        try:
            annot_df, abstract=fn.pubtator_annot(annot_df,PubID)
        except:
            abstract=fn.ret_abstract(PubID)
            if(abstract=='?'):
                    print("PMID=["+PubID+"] does not exist any more!")
                    continue # remove it from the output results in TRRUST
            else:
                    print("PubTator Response is not readable...!")
                    print("Try to annotate manually...")
        try:
            beCAS_lookup_full=fn.beCAS_lookup(PubID,query_id)
            beCAS_lookup=beCAS_lookup_full[['type','id','text','offset','length']]
            annot_df=pd.concat([annot_df,beCAS_lookup], ignore_index=True)
        except:
                print("beCAS Server error...!")

            #add surface similarity
        lookup_results=fn.lookup_annot(abstract,query_genes,query_id,lookup_ids)
        annot_df=pd.concat([annot_df,lookup_results],ignore_index=True)
        surface_annot=fn.surface_similarity(abstract, genes, query_genes, query_id,lookup_ids,single_id)
        annot_df=pd.concat([annot_df,surface_annot],ignore_index=True)
        annot_df=annot_df.drop_duplicates(subset=['id','offset'])

        annot_df=fn.multiple_taxonomy(annot_df, query_id)

        #depgraphs=fn.dep_parser('9000',abstract,annot_df,query_id,single_id,Positive,Negative)
        #node_topo=nx.topological_sort(depgraphs)
        #for node in depgraphs:
            #if depgraphs.out_degree(node)==0: #it's a leaf
              #  paths.append(nx.shortest_path(G, root, node))




        annot_df=annot_df.reset_index(drop=True)
        candidate_sentences, covered=fn.candidate_sentence(annot_df,abstract,query_id)

        if(len(candidate_sentences.index)==0):
                print('No co-existed sentences found in the abstract...!')
                continue
        target_sentences=[]
        for sentence in candidate_sentences.itertuples():
            obj=Rules(Positive,Negative,annot_df,covered,abstract,query_genes,query_id,sentence)
            depgraphs=fn.dep_parser('9000',sentence,annot_df,query_id,single_id,Positive,Negative,2)
            if(depgraphs):
                obj. multiplication_score(depgraphs, single_id)
            else:
                continue
            #obj.search_ranking()
            ranks.append(obj.rank)
            if(obj.rank!=0):
                target_sentences.append([sentence.sentence,obj.rank])

        rank_T1=sum(ranks)
        mode=''
        if(len(ranks)==0):
            continue
        if(rank_T1==0): sum_ranks.append(rank_T1)
        if(rank_T1>0): mode='positive'
        if(rank_T1<0): mode='negative'
        for sentence in target_sentences:
            print(str(single_id[0]) + "\t" + str(single_id[1]) + "\t" + single_gene[0]+ "\t" + single_gene[1] + "\t" + mode + "\t" + str(sentence[1]) + "\t" + sentence[0] + "\n")

        sum_ranks.append(rank_T1)


    rank_T2=sum(sum_ranks)
    if(len(sum_ranks)==0):
        print("There is no ranking value...!")
    if(rank_T2==0 and len(sum_ranks)!=0):
        print("The rank value is zero...!")
    if(rank_T2>0):
        print(str(single_id[0]) + "\t" + str(single_id[1]) + "\t" + single_gene[0]+ "\t" + single_gene[1] + "\t" + "positive" + "\t" + str(rank_T2) + "\n")
    if(rank_T2<0):
        print(str(single_id[0]) + "\t" + str(single_id[1]) + "\t" + single_gene[0]+ "\t" + single_gene[1] + "\t" + "negative" + "\t" + str(rank_T2) + "\n")
else:
    print("No PMID found for the interacion...!")













