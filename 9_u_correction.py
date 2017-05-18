#!/usr/bin/env python
'''
9_u_correction does the calculation for corrected <u> using "p", probability of having non-fixed states,
and protein alignments for each E.coli genes.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from os import listdir
from os.path import isfile, join
import math
import decimal as dec
import pandas as pd


def  main():
   
    # get list of alignments files
    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/alignments_nogaps"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    #try:
    
    # lists to create dataframe and csv
    id_list = []
    ulist = []
    nseq = []
    corr_list =[]
    p_list = []
    
    #read "p" from csv (p_calculation.py output)
    df = pd.read_csv('pia_calc.csv',sep="\t",tupleize_cols=1)
    idl = df['bid'].tolist()
    plist = df['pia'].tolist()
    
    #Go through alignment files and calculate corrected <u>
    for f in onlyfiles:
        
        if f == '.DS_Store':
            continue
        if (f.split('.')[0]) not in idl:
            continue
        print(f)
        
        #get p value from list
        pindex = idl.index(f.split('.')[0])
        p = plist[pindex]
        
        #calculate <u>
        u_nseq = calculate_u(f)
        #calculate corrected <u>
        corr_nseq = calculate_corr(f,p)
        umean= u_nseq[0]
        
        if umean != umean:
            continue
        
        #append values to lists
        nseq.append(u_nseq[1])
        ulist.append(umean)
        corr_list.append(corr_nseq)
        p_list.append(p)
        id_list.append(f.split('.')[0])
    
    #create dataframe
    df = pd.DataFrame({"b_id":id_list,"umean":ulist,"orth_count":nseq,"corr_u":corr_list,'p':p_list})

    #output to csv
    df.to_csv("u_correction.csv" , sep = "\t", index= False)
    
    
def calculate_u(filename):
    '''
    Calculate <u> for gene from protein sequence alignment file
    
    Arg: fasta file with protein alignments
    return: <u> and orth count as list
    
    '''
   
    records = list(SeqIO.parse(os.path.join('alignments_nogaps/',filename), "fasta"))
    
    # Go through records and save AA seq in seq_list
    seq_list = []
    
    for record in records:
        
        slist = list(str(record.seq))
        seq_list.append(slist)
    
    # Go through seq_list and make a list of lists with site AA
    site_list = []
    for y in range(0,len(seq_list[0])):
        #loop for site position (0 to protein length)
        site = []
        for x in range(0,len(seq_list)):
            #loop for all sequences
            site.append(seq_list[x][y])
        site_list.append(site)
            
    # Go through list of sites and calculate number of unique AA        
    u_list = []
    for site in site_list:
        
        aaset = set(site)
        if '-' in aaset:
            u_list.append(len(aaset)-1)
        else:
            u_list.append(len(set(site)))
     
    #return mean of unique AA count per sites   
    return([np.mean(u_list),len(seq_list)])


def calculate_corr(filename,p):
    '''
    Calculate corrected <u> 
    
    Args: protein alignment file, "p" probability to have non-fixed state
    return: corrected <u>
    
    '''
    
    #dec.getcontext().prec = 100
    records = list(SeqIO.parse(os.path.join('alignments_nogaps/',filename), "fasta"))
    
    # Go through records and save AA seq in seq_list
    seq_list = []
    
    for record in records:
        
        slist = list(str(record.seq))
        seq_list.append(slist)
    
    # Go through seq_list and make a list of lists with site AA
    site_list = []
    for y in range(0,len(seq_list[0])):
        #loop for site position (0 to protein length)
        site = []
        for x in range(0,len(seq_list)):
            #loop for all sequences
            site.append(seq_list[x][y])
        site_list.append(site)
            
    # Go through list of sites and calculate sigma for site      
    u_list = []
    sig_list = []
    
    #loop through sites in alignment
    for site in site_list:
        
        aaset = set(site)
        N=len(site)
        sigma = 0
        
        #loop through amino acid observed at site
        for aa in set(site):
            if aa == '-':
                continue
            #number of time the amino acid is observed "k"
            k = site.count(aa) 
            #probability that the amino acid is non-fixed "m" (poisson formula)
            m = ((dec.Decimal(dec.Decimal(p*N)**dec.Decimal(k))*dec.Decimal(np.exp(-(dec.Decimal(p*N)))))/dec.Decimal(math.factorial(k)))
            #probability that amino acid is fixed "r"
            r = 1-m
            #add probabilities to have fixed amino acid for every amino acid observe "sigma"
            sigma = sigma + r
        
        sig_list.append(sigma)

    #return mean of sigma for gene
    return(np.mean(sig_list))            
        
main()