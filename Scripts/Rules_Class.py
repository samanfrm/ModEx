#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from nltk.tokenize import word_tokenize
import numpy as np
import functions as fn

class Rules:
    def __init__(self,Positive,Negative,annot_df,covered,abstract,query_gene,query_id,sentence):
        self.Positive=Positive
        self.Negative=Negative
        self.annot_df=annot_df #list of all annotations for the abstract
        self.covered=covered #list of annotations covered by sentences including query genes
        self.abstract=abstract
        self.query_id=query_id
        self.query_gene=query_gene
        self.sentence=sentence
        self.rank=0
        self.tfrank=0
        self.tgrank=0
        self.direction=0 # +1 : TF - Target // -1 : Target - TF // 0 : not clear

    def search_ranking(self):
        word_tokens = word_tokenize(self.sentence.sentence)
        word_tokens = [w for w in word_tokens if (len(w)>1)]
        covered=self.sentence.covered
        query_list=self.annot_df.loc[covered,['id','text','offset']]
        ids=self.query_id[0]+self.query_id[1]
        query_list=query_list.loc[query_list["id"].isin(ids),:]
        target_list=query_list.loc[query_list["id"].isin(self.query_id[1]),:]
        voting_list=[]
        rank=0
        for target in target_list.itertuples():
            pos_event=0
            neg_event=0
            class_list=[]
            for item in word_tokens:
                try:
                    offset=self.sentence.start+self.sentence.sentence.index(item) #find the offset of item
                    if(item in self.Positive):
                        class_list.append([offset,1])
                        pos_event+=1
                    elif(item in self.Negative):
                        class_list.append([offset,-1])
                        neg_event+=1
                except:
                    continue
            if (pos_event==0 and neg_event==0):
                rank+=0 # when no positive/negative events found in sentence
                break
            if(pos_event==len(target_list.index) and neg_event==0):
                rank+=1 #when all events are found in positive events
                break
            if (neg_event==len(target_list.index) and pos_event==0):
                rank+=-1
                break #when all events are found in negative events

                #if non of above
                #find the closets event to target as the class
                #finally, we do voting to change the rank
            class_offset=np.asarray(class_list, dtype=np.float32)[:,0]
            dist=abs(class_offset-target.offset)
            voting_list.append(class_list[np.argmin(dist)][1])
            if(len(voting_list)==len(target_list.index)): # run when all of the targets in sentence assigned to a class -> voting
                if(sum(voting_list)>0):
                    rank+=1 #voting to positive
                    break
                elif(sum(voting_list)<0):
                    rank+=-1
                    break
                else:
                    rank+=0
                    break

        self.rank+=rank
        
    def multiplication_score(self, depgraphs, single_id):
        self.tfrank=fn.dep_score_source(depgraphs,single_id,self.Positive,self.Negative)
        self.tgrank=fn.dep_score_target(depgraphs,single_id,self.Positive,self.Negative)
        
        if(self.tgrank==0 and self.tfrank==0):
            self.rank+=0
            return
        if(self.tgrank==0 and self.tfrank!=0):
            if(self.tfrank>0):
                self.rank+=1
            else:
                self.rank+=-1
            return
        if(self.tgrank!=0 and self.tfrank==0):
            if(self.tgrank>0):
                self.rank+=1
            else:
                self.rank+=-1
            return

        if(self.tgrank!=0 and self.tfrank!=0):
            if(self.tgrank>0 and self.tfrank>0):
                self.rank+=1
                return

            if(self.tgrank<0 and self.tfrank<0):
                self.rank+=1
                return

            if(self.tgrank>0 and self.tfrank<0):
                self.rank+=-1
                return

            if(self.tgrank<0 and self.tfrank>0):
                self.rank+=-1
                return









