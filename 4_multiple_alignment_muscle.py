#!/usr/bin/env python
'''
4_multiple_alignment_muscle.py is used to do multiple alignment with MUSCLE on many 
fasta files (all E.coli genes). 

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
import re
from os import listdir
from os.path import isfile, join
import os


def main():
    
    mypath = "/RQusagers/egirard/Epistasis/orthologs_sub2"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    
    mypath2 = "/RQusagers/egirard/Epistasis/alignments_DNA"
    alignfiles = [f for f in listdir(mypath2) if isfile(join(mypath2, f))]
    
    #loop through fasta files to align
    for f in onlyfiles:
        if f in alignfiles:
            continue    
        if f == '.DS_Store':
            continue
        try:
            multiple_alignment(f)
        except:
            continue
        
def multiple_alignment (fasta_file):

    '''
    Output muscle command line for multiple alignment
    
    Arg: fasta file with sequences to align
    
    '''

    infile = os.path.join('orthologs_sub2/',fasta_file)
    outfile = os.path.join('alignments_DNA/',fasta_file)
    muscle_cline = 'muscle -in {} -out {} -maxiters 2'.format(infile,outfile)

    os.system(muscle_cline)
   


            
main()