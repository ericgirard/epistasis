#!/usr/bin/env python
'''
1_get_orthologs_seq.py create fasta files containing all orthologs sequences of E.coli 
genes within gammaproteobacteria. Those fastas can then be aligned and used for substitution 
rates calculation.

Author: Eric Girard

Serohijos Lab, Faculty of Medicine, Department of Biochemistry and Molecular Medicine, University de Montreal
'''

import pandas as pd
from bioservices.kegg import KEGG
#from Bio.KEGG import REST
import os.path
from os import listdir
from os.path import isfile, join


def main():
    
    #create csv with organism keggids 
    create_keggids_csv("ecoli_CDS.csv", "eco")
    
    #get list of genes files to know which one has already been downloaded (in case the script crash)
    mypath = "/Users/Eric/Dropbox/Stage_ete_2016/Epistasis/orthologs_fastas"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    
    #extract orthologs ids from list of genes keggids as dictionnary
    o_dict = extract_orthologs("ecoli_keggids.csv")

    #extract orthologs sequences and output to fasta for each genes
    extract_sequences(o_dict,onlyfiles)

  
def create_keggids_csv(filename, org):
    '''
    Extract keggids for an organism and save it to a csv file
        
        args: filename is the file containing gene name/ locus for all the organism genes
              org is the abrievation of the organism in kegg

    '''  
    
    #Open csv as panda dataframe (df)
    df = pd.read_csv(filename,sep="\t",tupleize_cols=1)
    gene_list = tuple(df['Locus'].tolist())
    bid_list = tuple(df['Locus tag'].tolist())
    kid_list = []
    
    k = KEGG()
    
    #find keggid for each genes
    for gene in bid_list:
        kstrg = (k.find(org, gene))
        kid_list.append(kstrg.split()[1])
    
    #create new df and save it to csv
    new_df = pd.DataFrame(columns=['gene','b_id','kegg_id'])
    new_df.gene = gene_list
    new_df.b_id = bid_list
    new_df.kegg_id = kid_list
   
    new_df.to_csv("ecoli_keggids.csv" , sep = "\t", index= False)
    
def extract_orthologs(filename):
    '''
    Create dictionnary with keggid as key and list of orthologs as value
        
        arg: csv with keggids
        return : dict with orthologs
    
    '''
    
    orthos_dict = {}
    k = KEGG()
    
    #get list of gammaproteobacteria from csv
    df = pd.read_csv(filename,sep="\t",tupleize_cols=1)
    df_gamma = pd.read_csv('gammaproteo.csv',sep="\t",tupleize_cols=1)
    gamma_list = df_gamma['KEGG'].tolist()
    
  
    #loop through keggid to get orthologs
    for keggid in df['kegg_id']:
        
        if keggid == "no":
            continue
            
        print(str(keggid))
        ortho_list = []
        
        #get orthologs on kegg
        data = k.get(keggid)
        dict_data = k.parse(data)
        
        if isinstance( dict_data, int ):
            continue
       
       #loop through kegg orthologs data and verify that organisms are gammaproteobacteria
        for key, value in dict_data['GENES'].items():
            
            if key.lower() in gamma_list:
               # print(key.lower(), value.split('(')[0].split())
                para_num =  len(value.split('(')[0].split())
                para_list = []
                
                for i in range (0,para_num):
                    #print(value.split('(')[0].split()[i])
                    para_list.append(key.lower() + ":" + value.split('(')[0].split()[i])
                
                ortho_list.append(para_list)
       
        orthos_dict[keggid] = ortho_list
        
    
    return orthos_dict
    
    
def extract_sequences(dict,flist):
    '''
    Get orthologs sequences on KEGG and write to a fasta file for each kegg id
        
        arg: dictionnary with keggid as key and orthologs as value (list)
        
    '''
    k = KEGG()
    
    ocount = {}
    
    #loop through orthologs dictionnary to get sequences from kegg
    for key, list in dict.items():
        #print(key)
        if (key + ".fas") in flist:
            print(key +" is already created !!!")
            continue
        
        #create string with sequences to write fasta file for each genes
        string = ""
        for x in range ( 0 , len(list)):
            for i in range ( 0 , len(list[x])):
                data_seq =  k.get(list[x][i], option = "ntseq", parse=True)
                string = string + data_seq + "\n"
                #print(data_seq)
        
        print("writing : " + key + ".fas")
        #write file
        with open(os.path.join('orthologs_fastas/',key + '.fas'), 'w') as f:
            read_data = f.write(string)
        f.closed
   
main()
