import os, sys, math, statistics, glob, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of descs and a list of seqs.
    """
    descs = []
    seqs = []
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespaces
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):  # Description line
                if sequence:
                    seqs.append(sequence)
                    sequence = ""
                descs.append(line[1:])
            else:
                sequence += line
        if sequence:
            seqs.append(sequence)
    return descs, seqs

def read_input(input_file):
    filepaths = []


    return filepaths

def link_genome_taxon(descs, taxonfile):

    return df

def taxon_summary()





input_pMSAs, 