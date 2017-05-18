#!/usr/bin/env python
'''
2_filter_seq_size.py removes orthologs sequences in fasta files being 15% longer or shorter 
than reference sequence ("eco").

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import os.path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os.path
from os import listdir
from os.path import isfile, join


def main():
    
    #get list of orthologs fastas file names
    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/orthologs_fastas"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    
    #loop through genes orthologs fasta file to filter sequences
    for f in onlyfiles:
        #try:
        if f == '.DS_Store':
            continue
        records = records_filter(f)
      
        write_fasta_file(records[1],records[0])
        
        #except:
            #print("nofile")

    
        
def records_filter(filename):
    '''
    Remove sequences 15% longer or shorter than reference sequence in fasta file
     
    Arg : fasta file name
    return : list of sequences as records and file name as a list
    
    '''
    
    seqlist =[]
    records = list(SeqIO.parse(os.path.join('orthologs_fastas/',filename), "fasta"))   
    refseq = ""
    fname = ""
    
    #look for reference sequence in file (ecoli)
    for record in records:
        if "eco:" in record.id:
            refseq = record.seq
            fname = record.id.split(":")[1]
    
    l = len(refseq)
    
    n = 0
    #loop through records and filter sequences by length
    for record in records:
        n = n + 1
        #print(record.id)
        #print(record.description)
        if (len(record.seq) != l) :
            if ((len(record.seq)/l)<=0.85) or ((len(record.seq)/l)>=1.15):
                #print(filename)
                #print((len(record.seq)/l))
                #print(n)
                continue
            else:
                seqlist.append(record)      
        else:
            seqlist.append(record)
            
    return [seqlist, fname]
        
def write_fasta_file(filename, records):
    '''
    Write records as fasta
    
    Args: file name to write, records to write
    '''
    print('Writing:{}'.format(filename+".fas"))
    output_handle = open(os.path.join('orthologs_sub2/',filename+".fas"), "w")
    SeqIO.write(records, output_handle, "fasta")
    output_handle.close()
        
        

main()
