#!/usr/bin/env python
'''
8_p_calculation.py calculate the probability "p" of observing a non-fixed state amino acid.
These "p" can further be used for correction of <u>.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

from bioservices.kegg import KEGG
import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from os import listdir
from os.path import isfile, join
import pandas as pd


def main():


    # get list of files
    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/alignments_nogaps/"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    
    df = pd.read_csv("pia_calc.csv",sep="\t",tupleize_cols=1)
    #lists to be added to panda df
    pia_all = df['pia'].tolist()
    bid_all = df['bid'].tolist()
    fract_all = df['fraction'].tolist()
    
    #loop through proteins alignments files(per gene)
    for f in onlyfiles:
        
        if (f.split('.')[0]) in bid_all:
            continue
        try: 
            if f == '.DS_Store':
                continue
            print(f)
  
            orgd = get_seq(f) #get dictionary of organisms and sequences
            pial = []
            total = 0 
            totcomp = 0
        
            #loop through organism dictionary for pairwise amino acid diversity calculation
            for key, value in orgd.items():
                total = total + len(value) # add count of sequences 
                if len(value) > 1:
                   # print(key)
                    totcomp = totcomp + len(value) # add count of sequences used in diversity calculation 
                    pia = calculate_pia(value) # diversity calculation
                    pial.append(pia)
                    '''
                    if pia > 0.1:
                        print(key)
                        print(pia)
                        for val in value:
                            print("".join(val))
                    '''
                
            bid_all.append(f.split('.')[0]) # add locus tag 
            fract_all.append(str((totcomp/total))) #add fraction of used sequences for the calculation
            pia_all.append(np.mean(pial)/6) #add average of pi (for species) divided by 6 to get p
            
            #create dataframe
            df = pd.DataFrame({'bid':bid_all,'pia':pia_all,'fraction':fract_all})
            #output dataframe to csv
            df.to_csv("pia_calc.csv" , sep = "\t", index= False)
            
        except Exception as e:
             print(str(e))
             continue
            
def get_seq(filename):
    '''
    Create dictionnary with species as keys and sequences as values for an alignment
    
    arg: filename with gene name
    return: organism dictionnary with sequences
    '''
    
    k = KEGG()
    records = list(SeqIO.parse(os.path.join('alignments_nogaps/',filename), "fasta"))
    
    idlist = []
    orglist = [] 
    seqlist = []
    orgdict = {}
    
    #go through sequences and search for organism name on kegg
    for record in records:
        
        idsplit = (record.id).split('_',1)
        id = idsplit[0] + ':' + idsplit[1]
        
        handle  = k.get(id)
        if isinstance( handle, int ):
            print(id)
            continue
            
        org = k.parse(handle)['ORGANISM']
        org = org.split()
        org = org[1] +" "+ org[2]
        seqlist.append(list(str(record.seq)))
        orglist.append(org)
        idlist.append(id)

    duplist = set(orglist)
    
    # create dict with organism as key and sequences for organism as values
    for org in duplist:
        indices = [i for i, x in enumerate(orglist) if x == org]
        seqs = []
        for e in indices:
            seqs.append(seqlist[e])
        orgdict[org] = seqs
        
    #print(orgdict)
    return orgdict
    
def calculate_pia(seq_list):
    '''
    Calculate pairwise amino acid diversity (pi) for all sequences in input (usually all from same specie)
    
    arg: list of sequences 
    return: mean of pairwise pi
    '''
    pia_list = []
    i=0
    #loop through specie sequences
    for seq in seq_list:
        i = i + 1 #add index to skip pair already compared
        #loop through specie sequences 
        for e in range(i,len(seq_list)):
            count = 0
            #loop through amino acid in sequences and do pairwise comparison for amino acids
            for x in range(0,len(seq)):
                if seq[x] != seq_list[e][x]:
                    count = count + 1 # add count if amino acid different
            pia_list.append(count/len(seq)) # normalize count by length of sequence 
    
    #return average across pairwise comparisons for specie
    return np.mean(pia_list)
    
main()