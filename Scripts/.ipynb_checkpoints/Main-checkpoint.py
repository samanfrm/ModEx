#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import sys
from Bio import Entrez
from Rules_Class import Rules
import functions as fn
from configCHIP import output_directory
from configCHIP import input_directory
from configCHIP import result_directory
#import configCHIP
#part=sys.argv[1]
part='part0'

Positive=[]
[Positive.append(line.strip().upper()) for line in open(input_directory+"Positive.txt")]
Negative=[]
[Negative.append(line.strip().upper()) for line in open(input_directory+"Negative.txt")]

genes_ents=input_directory + "ALL_Human_Genes_Info.csv" #NCBI
genes=pd.read_csv(genes_ents,sep=',',header=(0))
genes.fillna('', inplace=True)

lookup_ids=pd.read_csv(input_directory+"ncbi_id_lookup.csv",sep='\t',header=(0))

chip_ents=input_directory + "ChIPfilter-1K-950_ents.txt"
rels=pd.read_csv(output_directory+part+".csv",sep='\t',header=(0),dtype=object)
ents=pd.read_csv(chip_ents,sep='\t',header=(0))



CHIP_result=pd.DataFrame(columns=['srcuid','trguid','src_entrez','trg_entrez','srcname','trgname','find_pmid','all_pmids','mode','score','evi_pmid','evi_sent','report'])


for row in rels.iterrows():
        query_genes, query_id=fn.find_names(row,ents,genes)
        query_id=[int(query_id[0]),int(query_id[1])]
        print(query_genes)
        query_genes, query_id,single_gene,single_id =fn.make_query(genes,lookup_ids,query_genes,query_id)
        status=''
        try:
            myterm=fn.term_maker(single_gene,genes)
            ####  ESearch: Searching the Entrez databases
            Entrez.email="saman.farahmand001@umb.edu"
            handle=Entrez.esearch(db="pubmed", term=myterm, retmax=100000000)
            record=Entrez.read(handle)
        except:
            status+="Enterz Fetch Error|||"
            print(status)
            CHIP_result=CHIP_result.append({'srcuid':row[1]['srcuid'],'trguid':row[1]['trguid'],'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'find_pmid':None,'all_pmids':None,'mode':None,'score':None,'evi_pmid':None,'evi_sent':None,'report':status},ignore_index=True)
            continue
        PubIDs = record["IdList"]
        if(len(PubIDs)>0):
            sum_ranks=[]
            evi_pmids=[]
            evi_sentence=[]
            all_pmids=';'.join(PubIDs)
            for PubID in PubIDs:
                abstract=''
                ranks=[]
                annot_df=pd.DataFrame(columns=['type','id','text','offset','length'])
                try:
                    annot_df, abstract=fn.pubtator_annot(annot_df,PubID)
                except:
                    abstract=fn.ret_abstract(PubID)
                    if(abstract=='?'):
                        status+="PMID=["+PubID+"] does not exist any more|||"
                        continue # remove it from the output results in TRRUST
                    else:
                        status+="PMID=["+PubID+"] PubTator Response is not readable, Try to annotate manually|||"
                        #print(status)
#                try:
#                    beCAS_lookup_full=fn.beCAS_lookup(PubID,query_id)
#                    beCAS_lookup=beCAS_lookup_full[['type','id','text','offset','length']]
#                    annot_df=pd.concat([annot_df,beCAS_lookup], ignore_index=True)
#                except:
#                    status+="beCAS Server error|||"

                lookup_results=fn.lookup_annot(abstract,query_genes,query_id,lookup_ids)
                annot_df=annot_df.append(lookup_results)
#                surface_annot=fn.surface_similarity(abstract, genes, query_genes, query_id,lookup_ids,single_id)
#                annot_df=annot_df.append(surface_annot)
                annot_df=annot_df.drop_duplicates(subset=['id','offset'])

                annot_df=fn.multiple_taxonomy(annot_df, query_id)


                annot_df=annot_df.reset_index(drop=True)
                candidate_sentences, covered=fn.candidate_sentence(annot_df,abstract,query_id)
                if(len(candidate_sentences.index)==0):
                    status+="PMID=["+PubID+"] No co-existed sentences found in the abstract|||"
                    #print(status)
                    continue
                for sentence in candidate_sentences.itertuples():
                    obj=Rules(Positive,Negative,annot_df,covered,abstract,query_genes,query_id,sentence)
                    depgraphs=fn.dep_parser('9000',sentence,annot_df,query_id,single_id,Positive,Negative,2)
                    if(depgraphs):
                        try:
                            obj. multiplication_score(depgraphs, single_id)
                        except:
                            status+="PMID=["+PubID+"] dependency graph score error|||"
                    else:
                        status+="PMID=["+PubID+"] dependency graph co-occurance of single ids error|||"

                    #obj.search_ranking()
                    ranks.append(obj.rank)
                    if(obj.rank!=0):
                        evi_sentence.append('['+PubID+']'+sentence.sentence)
                        evi_pmids.append(PubID)
                if(len(ranks)!=0):
                    sum_ranks.append(sum(ranks))
            mode=''
            rank_T=sum(sum_ranks)
            if(rank_T>0):
                mode='positive'
            if(rank_T<0):
                mode='negative'

            evi_sentence=';'.join(evi_sentence)
            evi_pmids=';'.join(evi_pmids)
            CHIP_result=CHIP_result.append({'srcuid':row[1]['srcuid'],'trguid':row[1]['trguid'],'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'find_pmid':str(len(all_pmids)),'all_pmids':all_pmids,'mode':mode,'score':str(rank_T),'evi_pmid':evi_pmids,'evi_sent':evi_sentence,'report':status},ignore_index=True)
        else:
            status+="Not found any PMIDs for this interaction"
            print(status)
            CHIP_result=CHIP_result.append({'srcuid':row[1]['srcuid'],'trguid':row[1]['trguid'],'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'find_pmid':str('0'),'all_pmids':None,'mode':None,'score':None,'evi_pmid':None,'evi_sent':None,'report':status},ignore_index=True)


CHIP_result.to_csv(result_directory+part+"-CHIP-1K-950-result.csv",sep = '\t')