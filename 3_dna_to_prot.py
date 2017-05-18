#!/usr/bin/env python
'''
3_dna_to_prot.py is used to convert DNA sequences in fasta files to Amino acids (protein).

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from os import listdir
from os.path import isfile, join


def main():
        
        #get list of orthologs fastas file names
        mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/orthologs_sub2"
        onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
     
        for f in onlyfiles:
            if f == '.DS_Store':
                continue
            records = records_filter(f)
            write_fasta_file(f,records)
        
def records_filter(filename):
    '''
    Convert DNA sequence to protein
    
    Arg: fasta file name
    return: list of protein sequences as records
    
    '''
    
    seqlist =[]
    records = list(SeqIO.parse(os.path.join('orthologs_sub2/',filename), "fasta"))   
    
    n = 0
    #loop through records to convert DNA to proteins
    for record in records:
        n = n + 1
        #print(record.id)
        #print(record.description)
        newname = record.id.split()[0]
       # print(newname)
        record.description = ""
        record.id = newname
        prot = Seq.translate(record.seq,to_stop=True)
        record.seq = prot
        seqlist.append(record)
            
    return seqlist
        
def write_fasta_file(filename, records):
    '''
    Write records as fasta
    
    Args: file name to write, records to write
    '''
    print('Writing:{}'.format(filename))
    output_handle = open(os.path.join("orthologs_prot/",filename), "w")
    SeqIO.write(records, output_handle, "fasta")
    output_handle.close()
        
        

main()
