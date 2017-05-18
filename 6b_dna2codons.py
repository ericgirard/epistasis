#!/usr/bin/env python
'''
6b_dna3codons.py is used to do codons alignments from DNA sequences files and Protein alignment
for each E.coli genes.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
import os
import re
from os import listdir
from os.path import isfile, join


def main():

    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/orthologs_sub2"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    #try:
    
    for f in onlyfiles:    
        if f == '.DS_Store':
            continue
        if '.gb' in f:
            continue
        try:
            pal2nal(f)
            #gblock(f)
            #remove_spaces('{}.gb'.format(f))
        except:
            continue

def pal2nal (fasta_file):
    '''
    Create codons alignment with protein alignment and DNA sequences
    
    Arg: fasta file name to do codons alignment
    
    '''
    #input DNA fasta file
    dfile = os.path.join('orthologs_sub2/',fasta_file)
    #input protein alignment
    pfile = os.path.join('alignments/',fasta_file)
    #file to write codons alignment
    output =  os.path.join('alignments_codons/',fasta_file)
    
    pal2nal_cline = "perl /Users/Eric/Documents/pal2nal.v14/pal2nal.pl  {} {}  -output fasta >  {}".format(pfile,dfile,output)
    print(pal2nal_cline)
    
    os.system(str(pal2nal_cline))
    
        
def gblock(fasta_file):
    '''
    Output gblock command line to remove gaps from alignment
    
    Arg: fasta file with sequences alignment
    '''
    
    file = os.path.join('alignments_codons/',fasta_file)
    # clustalo_cline = ClustalOmegaCommandline(cl_path, infile=fasta_file +".fas",outfile=fasta_file+"_aln.fas",verbose=True, force=True)
    gblocks_cline = "/Users/Eric/Documents/Gblocks/Gblocks {} -b5=n -p=n -t=c -e=.gb".format(file)
    print(gblocks_cline)
    
    os.system(str(gblocks_cline))
    
  
def remove_spaces(file_name):
    '''
    Remove spaces from alignment processed by gblocks 
    
    Arg: fasta file with alignment to remove spaces
    
    '''
    
    list = []
    
    file_path = os.path.join('alignments_codons/',file_name)
    file_path2 = os.path.join('alignments_codons/',file_name)
    
    #Loop through file lines to remove spaces
    with open(file_path, 'r') as file:
        
        for line in file:
            s1 = line
            s2 = re.sub('[\s+]', '', s1)
            s2 = s2 + "\n" 
            list.append(s2)  
    
    #Write to file   
    file = open(file_path2, 'w')
    for line in list:
        file.write(line)
   
    file.close()

main()