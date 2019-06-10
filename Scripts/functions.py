#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import nltk
nltk.download('punkt')
import xml.etree.ElementTree as ET
import numpy as np
import re
from Bio import Entrez
from Bio import Medline
from nltk.tokenize import word_tokenize
from nltk.util import ngrams
import jellyfish
import becas
from nltk.parse.corenlp import CoreNLPDependencyParser
import networkx as nx
import requests
import string
import json
import sys
import os
from Rules_Class import Rules


def modex(query_id,parser_port,Positive,Negative,lookup_ids,genes,mesh='',email=''):
    
    CHIP_result=pd.DataFrame(columns=['src_entrez','trg_entrez','srcname','trgname','mode','score','evi_pmid','evi_sent'])
    try:
        query_genes=retrieve_symbol(query_id)
    except:
        print("No annotation has been retrieved associated with input Entrez IDs!")
        sys.exit(-1)

    query_genes, query_id,single_gene,single_id =make_query(genes,lookup_ids,query_genes,query_id)
    status=''
    if(mesh==''):
        mesh='humans'
    try:
        myterm=term_maker(single_gene,genes,mesh)
        ####  ESearch: Searching the Entrez databases
        Entrez.email=email
        handle=Entrez.esearch(db="pubmed", term=myterm, retmax=2000)
        record=Entrez.read(handle)
        PubIDs = record["IdList"]
    except:
        status+="Enterz Fetch Error|||"
        print("Enterz Fetch Error!")
        sys.exit(-1)
        #print(status)
        CHIP_result=CHIP_result.append({'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'mode':None,'score':None,'evi_pmid':None,'evi_sent':None},ignore_index=True)
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
                    annot_df, abstract=pubtator_annot(annot_df,PubID)
                except:
                    abstract=ret_abstract(PubID)
                    if(abstract=='?'):
                        status+="PMID=["+PubID+"] does not exist any more|||"
                        continue # remove it from the output results in TRRUST
                    else:
                        status+="PMID=["+PubID+"] PubTator Response is not readable, Try to annotate manually|||"
                        #print(status)
    #                try:
    #                    beCAS_lookup_full=beCAS_lookup(PubID,query_id)
    #                    beCAS_lookup=beCAS_lookup_full[['type','id','text','offset','length']]
    #                    annot_df=pd.concat([annot_df,beCAS_lookup], ignore_index=True)
    #                except:
    #                    status+="beCAS Server error|||"

                lookup_results=lookup_annot(abstract,query_genes,query_id,lookup_ids)
                annot_df=annot_df.append(lookup_results)
    #           surface_annot=surface_similarity(abstract, genes, query_genes, query_id,lookup_ids,single_id)
    #                annot_df=annot_df.append(surface_annot)
                annot_df=annot_df.drop_duplicates(subset=['id','offset'])

                annot_df=multiple_taxonomy(annot_df, query_id)

                annot_df=annot_df.reset_index(drop=True)
                candidate_sentences, covered=candidate_sentence(annot_df,abstract,query_id)
                if(len(candidate_sentences.index)==0):
                    status+="PMID=["+PubID+"] No co-existed sentences found in the abstract|||"
                    #print(status)
                    continue
                for sentence in candidate_sentences.itertuples():
                    obj=Rules(Positive,Negative,annot_df,covered,abstract,query_genes,query_id,sentence)
                    depgraphs=dep_parser(parser_port,sentence,annot_df,query_id,single_id,Positive,Negative,2)
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
            evi_sentence='|||'.join(evi_sentence)
            evi_pmids=';'.join(evi_pmids)
            CHIP_result=CHIP_result.append({'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'mode':mode,'score':str(rank_T),'evi_pmid':evi_pmids,'evi_sent':evi_sentence},ignore_index=True)
            print(CHIP_result["srcname"].values,',',CHIP_result["trgname"].values)
            print(CHIP_result["mode"].values)
            print(CHIP_result["score"].values)
            print(CHIP_result["evi_pmid"].values)
    else:
        status+="Not found any PMIDs for this interaction"
        #print(status)
        CHIP_result=CHIP_result.append({'src_entrez':single_id[0],'trg_entrez':single_id[1],'srcname':single_gene[0],'trgname':single_gene[1],'mode':None,'score':None,'evi_pmid':None,'evi_sent':None},ignore_index=True)

    return CHIP_result






def find_names(row,ents,genes):
    query_genes=[]
    query_id=[]
    query_genes=ents[ents['uid'].isin([row[1]['srcuid'],row[1]['trguid']])]['name'].tolist()
    query_id=ents[ents['uid'].isin([row[1]['srcuid'],row[1]['trguid']])]['id'].tolist()
    for i in range(2):
        if(np.isnan(query_id[i])): #if not found in ChipAtlas then query to NCBI
            query_id[i]=genes[genes['Symbol']==query_genes[i]]['GeneID'].values[0]
    query_id =[str(int(x)) for x in query_id]
    return query_genes, query_id


def process_text(text):
    reach_text_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/text'
    data = {'text': text.encode('utf-8')}
    try:
            res = requests.post(reach_text_url, data)
    except:
            return
    json_str = res.content
    json_str = json_str.decode('utf-8')
    json_str = json_str.replace('frame-id', 'frame_id')
    json_str = json_str.replace('argument-label', 'argument_label')
    json_str = json_str.replace('object-meta', 'object_meta')
    json_str = json_str.replace('doc-id', 'doc_id')
    json_str = json_str.replace('is-hypothesis', 'is_hypothesis')
    json_str = json_str.replace('is-negated', 'is_negated')
    json_str = json_str.replace('is-direct', 'is_direct')
    json_str = json_str.replace('found-by', 'found_by')
    try:
        json_dict = json.loads(json_str)
    except:
        print('json load error')
    return json_dict
def REACH_apiv2(candidate_sentence,query_genes,query_id,Positive,Negative):
    for sent in candidate_sentence.itertuples():
        jason_dic=process_text(sent.sentence)
        flag=0
        for frame in jason_dic['events']['frames']:
            if(frame['type'] in ['regulation','activation']):
                if(frame['subtype'] in ['negative-regulation','negative-activation']):
                    return flag, 'negative'
                if(frame['subtype'] in ['positive-regulation','positive-actication']):
                    return flag, 'positive'
    return flag, None

def REACH_api(candidate_sentence,query_genes,query_id,Positive,Negative):
    for sent in candidate_sentence.itertuples():
        jason_dic=process_text(sent.sentence)
        tf_hit=0
        trg_hit=0
        flag=0
        entities={}
        for entity in jason_dic['entities']['frames']:
            id_list=[]
            for xref in entity['xrefs']:
                id_list.append(xref['id'])
            entities[entity['frame_id']]=id_list

        for frame in jason_dic['events']['frames']:
                    entrez_id=[]
                    for argument in frame['arguments']:
                        if(argument['argument-type']=='entity'):
                            ids=entities[argument['arg']]
                            entrez_id=convert_uni_2_entrez_list(ids,query_id)

                        if(argument['argument-type']=='event'):
                            for quer in query_genes[0]:
                                if(argument['text'].find(quer)!=-1):
                                    entrez_id.append(query_id[0])
                            for quer in query_genes[1]:
                                if(argument['text'].find(quer)!=-1):
                                    entrez_id.append(query_id[1])


                        if(any(map(lambda word: word in entrez_id ,query_id[0])) and any(map(lambda word: word in entrez_id ,query_id[1]))):
                            trg_hit=1
                            tf_hit=1
                    if(trg_hit==1 and tf_hit==1):
                        flag=1
                        if(frame['type']=='regulation'):
                            if(frame['subtype']=='negative-regulation'):
                                return flag, 'negative'
                            if(frame['subtype']=='positive-regulation'):
                                return flag, 'positive'
    return flag, None

def send_query(base_url,text, service_endpoint='drum', query_args=None):
    if service_endpoint in ['drum', 'drum-dev', 'cwms', 'cwmsreader']:
        url = base_url + service_endpoint
    if query_args is None:
        query_args = {}
    query_args.update({'input': text})
    res = requests.post(url, query_args)
    return res.text


def get_xml(html, content_tag='ekb', fail_if_empty=False):
    cont = re.findall(r'<%(tag)s(.*?)>(.*?)</%(tag)s>' % {'tag': content_tag},
                      html, re.MULTILINE | re.DOTALL)
    if cont:
        events_terms = ''.join([l.strip() for l in cont[0][1].splitlines()])
        if 'xmlns' in cont[0][0]:
            meta = ' '.join([l.strip() for l in cont[0][0].splitlines()])
        else:
            meta = ''
    else:
        events_terms = ''
        meta = ''

    if fail_if_empty:
        assert events_terms != '',\
            "Got empty string for events content from html:\n%s" % html

    header = ('<?xml version="1.0" encoding="utf-8" standalone="yes"?><%s%s>'
              % (content_tag, meta))
    footer = '</%s>' % content_tag
    return header + events_terms.replace('\n', '') + footer


def TRIPS_apiv2(candidate_sentence,query_genes,Positive,Negative):
    for sent in candidate_sentence.itertuples():
        URL='http://trips.ihmc.us/parser/cgi/'
        html = send_query(URL,sent.sentence,'drum')
        xml =  get_xml(html)
        tree = ET.ElementTree(ET.fromstring(xml))
        events=tree.findall('EVENT')
        flag=0
        for e in events:
           event_type = e.find('type').text
           mode=event_type.replace('ONT::','')
           arg1=e.find('arg1')
           if(arg1 is None):
               continue
           arg1=arg1.find('text').text
           arg2=e.find('arg2')
           if(arg2 is None):
               continue
           arg2=arg2.find('text').text
           flag=1
           if(mode in Positive):
               return flag, 'positive'
           if(mode in Negative):
               return flag, 'negative'

    return flag, None

def TRIPS_api(candidate_sentence,query_genes,Positive,Negative):
    for sent in candidate_sentence.itertuples():
        URL='http://trips.ihmc.us/parser/cgi/'
        html = send_query(URL,sent.sentence,'drum')
        xml =  get_xml(html)
        tree = ET.ElementTree(ET.fromstring(xml))
        events=tree.findall('EVENT')
        flag=0
        for e in events:
           event_type = e.find('type').text
           mode=event_type.replace('ONT::','')
           arg1=e.find('arg1')
           if(arg1 is None):
               continue
           arg1=arg1.find('text').text
           for quer in query_genes[0]:
               if(arg1.find(quer)!=-1):
                   arg2=e.find('arg2')
                   if(arg2 is None):
                       continue
                   arg2=arg2.find('text').text
                   for quer2 in query_genes[1]:
                       if(arg2.find(quer2)!=-1):
                          flag=1
                          if(mode in Positive):
                              return flag, 'positive'
                          if(mode in Negative):
                              return flag, 'negative'

    return flag, None

def REACH_extraction(candidate_sentence,query_genes,Positive,Negative):
    all_statements=[]
    for sent in candidate_sentence.itertuples():
        reach_processor = reach.process_text(sent.sentence)
        if reach_processor is not None:
            all_statements += reach_processor.statements

        for st in all_statements:
            st=str(st)
            event1=st[st.find("(")+1:st.rfind(")")]
            src=event1.split(',')[0].replace('()','').strip().upper()
            trg=event1.split(',')[1].replace('()','').strip().upper()
            mode=st.replace('('+event1+')','').strip().upper()
            if((src in query_genes[0]) and trg in query_genes[1]):
               max_pos=0
               max_neg=0
               for term in Positive:
                   score=jellyfish.jaro_distance(term,mode)
                   if(score>max_pos):
                       max_pos=score
               for term in Negative:
                   score=jellyfish.jaro_distance(term,mode)
                   if(score>max_neg):
                       max_neg=score
               if max_pos>max_neg:
                    return 'positive'
               else:
                    return 'negative'




def multiple_taxonomy(annot_df, query_id):

    unique_list=list(set(annot_df['text']))
    for gene in unique_list:
        corresponding=annot_df.loc[annot_df['text']==gene,'id'].tolist()
        unique_ids=list(set(corresponding))
        if(len(unique_ids)>1):
            for s_id in unique_ids:
                if(s_id in query_id[0]) or (s_id in query_id[1]):
                    continue
                else:
                    annot_df=annot_df[annot_df['id']!=s_id]
    return annot_df




def dep_score_source(depgraphs,single_id,positive_list,negative_list):
    score=0
    source=str(single_id[0])
    events=[]
    ancs=nx.ancestors(depgraphs,source)
    for node in ancs:
        if(node in positive_list):
            events.append([node,1])
            continue
        if(node in negative_list):
            events.append([node,-1])
            continue
    if(len(events)>0):
        lens=[]
        for event in events:
            event_node=event[0] #name of node
            lens.append(nx.shortest_path_length(depgraphs,event_node,source))
        ind=np.argmin(lens)
        if(events[ind][1]==1):
            score=1
        else:
            score=-1

    return score

def dep_score_target(depgraphs,single_id,positive_list,negative_list):
    score=0
    target=str(single_id[1])
    events=[]
    ancs=nx.ancestors(depgraphs,target)
    for node in ancs:
        if(node in positive_list):
            events.append([node,1])
            continue
        if(node in negative_list):
            events.append([node,-1])
            continue
    if(len(events)>0):
        lens=[]
        for event in events:
            event_node=event[0] #name of node
            lens.append(nx.shortest_path_length(depgraphs,event_node,target))
        ind=np.argmin(lens)
        if(events[ind][1]==1):
            score=1
        else:
            score=-1

    return score



def dep_parser(port,abstract,annot_df,query_id,single_ids,positive_list,negative_list,code):

    dep_parser=CoreNLPDependencyParser(url='http://localhost:'+port)
    if(code==1): #abstract
        sentences=nltk.sent_tokenize(abstract)
        offset=0
    else: #single sentence=2
        sentences=[abstract.sentence]
        offset=abstract.start
        single_ids=[str(x) for x in single_ids]
    dep_graph=nx.DiGraph()
    for sentence in sentences:
        end=offset+len(sentence)
        if(code==1):
            covered=annot_df.loc[(annot_df['offset']>=offset) & (annot_df['offset']<=end),:]
            ids=covered['id'].tolist()
        else:
            covered=abstract.covered
            covered=annot_df.loc[covered,:]
            covered=covered.drop_duplicates(subset=['id','offset'])
            covered=multiple_taxonomy(covered, query_id)
            ids=covered['id'].tolist()

        #normalize covered df by mapping query ids to single_ids
        #make two lists 1) query genes 2) non-query genes
        covered=covered.reset_index(drop=True)
        covered_query=pd.DataFrame(columns=['type','id','text','offset','length'])
        for i in covered.index.tolist():
            if(covered.loc[i].id in query_id[0]):
                covered_query=covered_query.append({'type':'Gene','id':str(single_ids[0]),'text':covered.loc[i].text.upper(),'offset':covered.loc[i].offset,'length':covered.loc[i].length},ignore_index=True)
                continue
            if(covered.loc[i].id in query_id[1]):
                covered_query=covered_query.append({'type':'Gene','id':str(single_ids[1]),'text':covered.loc[i].text.upper(),'offset':covered.loc[i].offset,'length':covered.loc[i].length},ignore_index=True)
                continue
        covered_nonquery=covered[(covered['id']!=single_ids[0]) & (covered['id']!=single_ids[1])]


        #first remove query genes
        #then remove non-query
        if (any(map(lambda word: word in ids ,query_id[0])) or any(map(lambda word: word in ids,query_id[1]))):
            for covered_gene in covered_query.itertuples():
                gene_id=covered_gene.id
                gene_text=covered_gene.text.upper()
                gene_text=re.escape(gene_text)
                if(re.search(r'\b('+gene_text+r')\b',sentence)):
                    sentence=re.sub(r'\b('+gene_text+r')\b',gene_id, sentence)
                    continue

            for covered_gene in covered_nonquery.itertuples():
                gene_id=covered_gene.id
                gene_text=covered_gene.text.upper()
                gene_text=re.escape(gene_text)
                if(re.search(r'\b('+gene_text+r')\b',sentence)):
                    sentence=re.sub(r'\b('+gene_text+r')\b',gene_id, sentence)
                    continue

            #remove punctuations
            tokens=word_tokenize(sentence)
#            new_tokens=[]
#            for tok in tokens:
#                if(len(tok)!=1):
#                    new_tokens.append(tok)

            #word_tokenize
            #replace string with substring of either src/target with the ids
#            tokens=new_tokens
            for query in single_ids:
                for i in range(0,len(tokens)):
                    if (tokens[i].find(query)!=-1 and len(tokens[i])!=len(query) and (tokens[i].isdigit()!=True)):
                        tokens[i]=query
            sentence=' '.join(tokens)

            if( ((not(single_ids[0] in tokens)) or (not(single_ids[1] in tokens))) and code==2):
                return '' # when code=2 and both single query are not in sentence


            dep_result=dep_parser.raw_parse(sentence)
            dep=dep_result.__next__()
            dep_list=list(dep.triples())
            for dep_rel in dep_list:
                src=[]
                trg=[]
                weig=1
                if( (dep_rel[0][0] in string.punctuation) or (dep_rel[2][1] in string.punctuation) ): # a node is punctuation
                    continue
                for i in range(0,3,2):
                    name=dep_rel[i][0]
                    mode=''
                    is_gene='0'
                    if(name in positive_list):
                        mode='positive'
                        weig+=1
                    if(name in negative_list):
                        mode='negative'
                        weig+=1
                    if(name in ids):
                        is_gene='1'
                        weig+=1
                    if(i==0):
                        src=[name,mode,is_gene]
                    else:
                        trg=[name,mode,is_gene]
                    dep_graph.add_node(name,label=name,sign=mode,isgene=is_gene)
                if(not dep_graph.has_edge(src[0],trg[0])):
                    dep_graph.add_edge(src[0],trg[0],weight=weig)
        offset=end+1

    return dep_graph


def term_maker(query_genes,genes,mesh='humans'):
    myterm=mesh+'[MeSH Terms]'
    for j in range(0,len(query_genes)):
        myterm+=' AND '
        myterm+='('
        flag=genes.loc[genes["Symbol"]==query_genes[j],"Synonyms"].tolist()
        syn=[]
        if(len(flag)>0 and flag[0]!='-') :
            syn=flag[0].split('|')
        syn.append(query_genes[j])
        for i in range(0,len(syn)):
            myterm+=syn[i]+'[sym]'
            if(i!=(len(syn)-1)):
                myterm+=' OR '
        myterm+=')'
    return myterm


def make_query(genes,lookup_ids,query_genes,query_id):
    genes_final=[[],[]]
    single_id=query_id
    ids_final=[[],[]]
    for i in range(2):
        flag=genes.loc[genes["Symbol"]==query_genes[i],"Synonyms"].tolist()
        syn=[]
        if(len(flag)>0 and flag[0]!='-'):
            syn=flag[0].split('|')
        syn.append(query_genes[i])
        genes_final[i]=syn
        if(np.isnan(query_id[i])): #if not found in ChipAtlas then query to NCBI
            query_id[i]=genes[genes['Symbol']==query_genes[i]]['GeneID'].values[0]
        alt_id=[]
        for j in syn: # extract ids for each synonym
            alt_id=alt_id+lookup_ids.loc[lookup_ids["symbol"]==j,"associated_ids"].tolist()[0].split('|')
        ids_final[i]=list(set(alt_id))

    return genes_final, ids_final, query_genes, single_id


def lookup_annot(abstract,query_genes,query_id,lookup_ids):
    # add lookup table results here to annot_df
    # send abstract to function
    # return results and append to annot_df
    lookup_annot=pd.DataFrame(columns=['type','id','text','offset','length'])
    syms=query_genes[0]+query_genes[1]
    symbols=syms
    #symbols=[]
    for symbol in symbols:
        for occur in re.finditer(r'\b('+symbol.upper()+r')\b',abstract):
            start=occur.start()
            length=len(symbol)
            ids=lookup_ids.loc[lookup_ids["symbol"]==symbol,"associated_ids"].tolist()
            if(len(ids)==0):
                return
            ids=ids[0].split('|')
            for _id in ids:
                lookup_annot=lookup_annot.append({'type':'Gene','id':str(_id),'text':symbol.upper(),'offset':start,'length':length},ignore_index=True)
    #remove overlaps based on ID and offset
    return lookup_annot


def pubtator_annot(annot_df,PubID):
    url = "https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/Gene/"+PubID+"/BioC"
    myResponse = requests.get(url)
    abstract=''
    if(myResponse.ok):
        POS_xml=myResponse.content.decode("utf-8-sig")
        #parser = ET.XMLParser(encoding="utf-8")
        tree = ET.ElementTree(ET.fromstring(POS_xml))
        flag=False
        for annotation in tree.findall('.//document/passage/annotation'):
            flag=True
            infon_type = annotation.findall('infon')[0].text
            infon_id = annotation.findall('infon')[1].text.split(';')
            for ID in infon_id:
                text = annotation.find('text').text.upper()
                offset = int(annotation.find('location').attrib['offset'])
                length = int(annotation.find('location').attrib['length'])
                annot_df=annot_df.append({'type':infon_type,'id':ID,'text':text,'offset':offset,'length':length},ignore_index=True)
                abstract=''
                if(flag):
                    for passage in tree.findall('.//document/passage'):
                        abstract+=passage.find('text').text.upper()+' '
    return annot_df, abstract

#find co-occurance of TF and Target members
def candidate_sentence(annot_df,abstract,query_id):
    sent_text = nltk.sent_tokenize(abstract)
    offset=0
    covered=pd.DataFrame()
    candidate_sentences=pd.DataFrame(columns=['sentence','start','end','covered']) # our target sentences with list covered gene annotations including query genes
    for sentence in sent_text:
        end=offset+len(sentence)
        covered=annot_df.loc[(annot_df['offset']>=offset) & (annot_df['offset']<=end),:]
        ids=covered['id'].tolist()
        if(any(map(lambda word: word in ids ,query_id[0])) and any(map(lambda word: word in ids,query_id[1]))): # query_id will be a list of ids not just 2
            candidate_sentences=candidate_sentences.append({'sentence':sentence,'start':offset,'end':end,'covered':covered.index.tolist()},ignore_index=True)
        offset=end+1
    return candidate_sentences, covered

def ret_abstract(pmid):
    Entrez.email="saman.farahmand001@umb.edu"
    Entrez.api_key="a41996d7fa7d12ce66f7ae12735b97b43708"
    handle = Entrez.efetch(db="pubmed", id=int(pmid), rettype="medline", retmode="text",) # See medline format table
    records = Medline.parse(handle)
    record = list(records)[0]
    abst=record.get("AB","?")
    return abst.upper()

def get_ngrams(text, n ):
    n_grams = ngrams(word_tokenize(text), n)
    return [ ' '.join(grams) for grams in n_grams]


def surface_similarity(abstract, genes, query_genes, query_ids,lookup_ids,single_id):
    candidates=[[],[]]#list names TF
    lookup_annot=pd.DataFrame(columns=['type','id','text','offset','length'])
    sentences=nltk.sent_tokenize(abstract)
    for i in range(0,len(query_genes)):
        candidates[i]=candidates[i]+query_genes[i]
        for symbol in query_genes[i]:
            ids=lookup_ids.loc[lookup_ids["symbol"]==symbol,"associated_ids"].tolist()
            if(len(ids)==0):
                return
            ids=ids[0].split('|')
            for _id in ids:
                desc=genes.loc[genes["GeneID"]==int(_id),"description"].tolist()[0].upper()
                candidates[i].append(desc)
        candidates[i]=list(set(candidates[i]))
        offset=0
        for sent in sentences:
            nlen=range(1,len(sent.split())+1)
            names=[]
            scores=[]
            for j in nlen:
                grams = get_ngrams(sent, j)
                for gram in grams:
                    words=candidates[i]
                    for word in words:
                        score=jellyfish.jaro_distance(word,gram)
                        names.append(gram)
                        scores.append(score)
                        #score=distance.get_jaro_distance(word, grams, winkler=True, scaling=0.1)
            ziip=list(zip(names, scores))
            max_entry = ('', 0)
            for entry in ziip:
                if entry[1] > max_entry[1]:
                    max_entry = entry
            if(max_entry[1]>= 0.9):
                maxent=re.escape(max_entry[0])
                #maxent=max_entry[0]
                try:
                    for occur in re.finditer(r'\b('+maxent+r')\b',sent):
                        start=offset+occur.start()
                        lookup_annot=lookup_annot.append({'type':'Gene','id':str(single_id[i]),'text':maxent.upper(),'offset':start,'length':len(maxent)},ignore_index=True)
                except:
                    continue
            offset+=len(sent)+1

    return lookup_annot

def beCAS_lookup(pmid, uidlist):
    becas.email = 'samanfm@gmail.com'
    becas.SEMANTIC_GROUPS = ('PRGE')
    results = becas.annotate_publication(int(pmid),groups={'PRGE':True})
    lookup_annot=pd.DataFrame(columns=['type','id','text','offset','length','description'])
    entities=[[],[]]
    entities[0]=results['entities_title']
    entities[1]=results['entities_abstract']
    gap=0
    for i in range(0,2): # one for title and one for abstract
        _entities=entities[i]
        gap=0
        if(i==1): gap=len(results['title'])
        for entity in _entities:
            chunks=entity.split('|')
            text=chunks[0]
            uniprotID=chunks[1].split(';')
            for _uniID in uniprotID:
                try:
                    segments=_uniID.split(':')
                    #desc=ids[_uniID]['name']
                    query_id=segments[1]
                    offset=gap+int(chunks[2])
                    EntrezID = convert_uni_2_entrez(query_id, uidlist)
                    if(any(map(lambda word: word in EntrezID ,uidlist[0])) and any(map(lambda word: word in EntrezID ,uidlist[1]))):
                        continue
                    for _id in EntrezID:
                        try:
                            lookup_annot=lookup_annot.append({'type':'Gene','id':str(_id),'text':text.upper(),'offset':offset,'length':len(text)},ignore_index=True)
                        except:
                            continue
                except:
                    continue

    return lookup_annot

def convert_uni_2_entrez(uniprot,uidlist):
    """Convert Uniprot Id to Entrez Id"""
    idlist=uidlist[0]+uidlist[1]
    trglist=[]
    Entrez.email="saman.farahmand001@umb.edu"
    convert=Entrez.esearch(db="gene", term=uniprot, retmax=10000)
    convert_records=Entrez.read(convert)
    geneId=list(set(convert_records["IdList"]))
        # check to see if more than one result is returned
        # if you have more than more result then check which Entrez Id returns the same uniprot Id entered.
    for _id in geneId:
        if(_id in idlist):
            trglist.append(_id)

    return trglist


def convert_uni_2_entrez_list(uniprotlst,uidlist):
    """Convert Uniprot Id to Entrez Id"""
    idlist=uidlist[0]+uidlist[1]
    trglist=[]
    Entrez.email="saman.farahmand001@umb.edu"
    geneId=[]
    for uniprot in uniprotlst:
        convert=Entrez.esearch(db="gene", term=uniprot, retmax=10000)
        convert_records=Entrez.read(convert)
        geneId+=list(set(convert_records["IdList"]))
        # check to see if more than one result is returned
        # if you have more than more result then check which Entrez Id returns the same uniprot Id entered.
    for _id in geneId:
        if(_id in idlist):
            trglist.append(_id)

    return trglist


def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene",id=",".join(str(x) for x in id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        sys.exit(-1)
    Entrez.email="samanfm@gmail.com"
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)

    return annotations

def retrieve_symbol(id_list):
    annotation=retrieve_annotation(id_list)
    gene_symbol=[]
    for gene_data in annotation["DocumentSummarySet"]["DocumentSummary"]:
        #gene_id = gene_data["Id"]
        gene_symbol.append(str(gene_data["NomenclatureSymbol"]))
    return gene_symbol
