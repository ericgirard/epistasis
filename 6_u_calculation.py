#!/usr/bin/env python
'''
6_u_calculation.py does the calculation for <u>, the average mutational usage, for each 
E.coli genes.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from os import listdir
from os.path import isfile, join
import pandas as pd


def  main():

    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/alignments_DNA2_nogaps"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    #try:
    
    #lists of bid, <u> and orthologs count
    id_list = []
    ulist = []
    nseq = []
    
    #loop through alignments files and calculate u and orthologs count
    for f in onlyfiles:
        print(f)
        if f == '.DS_Store':
            continue
        # calculate <u>
        u_nseq = calculate_u(f)
        umean= u_nseq[0]
        
        if umean != umean:
            continue
        # add orth count to list
        nseq.append(u_nseq[1])
        # add <u> to list
        ulist.append(umean)
        # add bid to list
        id_list.append(f.split('.')[0])
        dndsl.append(dnds_gene)
        break
    
    #create dataframe
    df = pd.DataFrame({"b_id":id_list,"umean":ulist,"orth_count":nseq})

    #output dataframe to csv
    df.to_csv("u_values.csv" , sep = "\t", index= False)
    
    
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
        
main()