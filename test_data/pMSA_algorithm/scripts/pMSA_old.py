import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# original pMSA algorithm that does all vs all symmetric pairing based on file name
def read_orthologues(files):
    #pass the list of file paths of orthologue sequences that have been aligned 
    MSA = pd.DataFrame()
    #fill dataframe with sequences and ID's
    for file in files:
        #split the gene name only from the file path, need to double check
        orthologue = file.split('/')[-1].split('.')[0]
        for record in SeqIO.parse(file, "fasta"):
            #create an entry with index species and column gene name. Append AA length ':len'
            MSA.loc[str(record.description), orthologue] = str(record.seq)
            #need to deal with NaN strings. "df.fillna('', inplace=True)"
            #maybe have unique string identifier if NaN to deal with later cleanup. hard to know which ones 
    return MSA

def pairedMSA(MSA):
#function that inputs the MSA dataframe and outputs a large .fasta of all pMSA 
#columns of MSA are orthologue names, rows are isolates
    pMSAs = pd.DataFrame()
    pairs = {} #dictionary that holds keys of added combinations of orthologues
    for orthologue_1 in MSA.columns:
        for orthologue_2 in MSA.columns:
            #make sure this combination of homolgoues has not been created before
            if (orthologue_1 != orthologue_2) and (not orthologue_1+orthologue_2 in pairs) and (not orthologue_2+orthologue_1 in pairs):
                pairs[orthologue_1+orthologue_2] = 1 #add new pair to dictionary set value to 1
                #calculate first chain length. Finds first string and uses this for length
                chain_length = []
                for index in MSA.index:
                    if type(MSA.loc[index, orthologue_1]) == str: #assure not NaN
                        chain_length = str(len(MSA.loc[index, orthologue_1]))
                        break
                #Column name orthologue_1__orthologue_2-AA# where # is the AA length of first chain
                pMSAs[orthologue_1+'__'+orthologue_2+'--AA'+chain_length] \
                    = MSA[orthologue_1]+MSA[orthologue_2]
    return pMSAs

def pMSAs_to_fasta(pMSAs):
    #this is going to input a dataframe with pMSAs and output a .fasta file
    #to the output directory for each A and B orthologue pair
    #ignores NaN entries
    for column in pMSAs.columns:
        #output file name is orthologue_1__orthologue_2-AA#.fasta where # is the AA length of first chain
        ofile = open(column+'.fasta', 'w')
        for index in pMSAs.index:
            if type(pMSAs.loc[index, column]) == str: #ignore NaN entries 
                #each .fasta entry has the format description orthologue_1__orthologue_2-AA#
                ofile.write('>'+index+' '+column+'\n'+pMSAs.loc[index,column]+'\n')
        ofile.close()
    return 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#input comma separated list of .fasta file paths as strings in Bash
MSA = read_orthologues(sys.argv[1:])
MSA.to_csv("run_MSA.csv")
pMSAs = pairedMSA(MSA) #generate df of pMSAs
pMSAs.to_csv("run_pMSA.csv")
pMSAs_to_fasta(pMSAs) #generate .fasta output



