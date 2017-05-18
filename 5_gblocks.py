#!/usr/bin/env python
'''
5_gblocks.py is used to remove gaps from alignments with Gblock and also remove
spaces in sequences processed by Gblock.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
import os
import re
from os import listdir
from os.path import isfile, join


def main():
    
    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/alignments_DNA2"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    #try:
    
    for f in onlyfiles:    
        if f == '.DS_Store':
            continue
        if '.gb' in f:
            continue
        gblock(f)
        remove_spaces('{}.gb'.format(f))

   # gblock('output.fas')
   # remove_spaces('output.fas.gb')



def gblock(fasta_file):
    '''
    Output gblock command line to remove gaps from alignment
    
    Arg: fasta file with sequences alignment
    '''
    
    file = os.path.join('alignments_DNA2/',fasta_file)
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
    
    file_path = os.path.join('alignments_DNA2/',file_name)
    file_path2 = os.path.join('alignments_DNA2_nogaps/',file_name)
    
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