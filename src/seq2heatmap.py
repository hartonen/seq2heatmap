#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import argparse
import csv
import numpy as np

from Bio import SeqIO

from matplotlib import pyplot as plt
from matplotlib import colors

def seq2heatmap():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()

    #MANDATORY PARAMETERS
    parser.add_argument("reads",help="Full path to the file containing the input reads in fasta-format.",type=str)
    parser.add_argument("outfile",help="Full path to the output file (png).",type=str)

    #OPTIONAL PARAMETERS
    parser.add_argument("--N",help="Number of first reads read in (default=1000).",type=int,default=1000)
    parser.add_argument("--L",help="Length of sequences (default=170).",type=int,default=170)
    parser.add_argument("--title",help="Plot title (default=None).",type=str,default=None)
    
    args = parser.parse_args()

    cols = {1:'green',2:'blue',3:'yellow',4:'red'} #A=1 C=2 G=3 T=4

    cvr = colors.ColorConverter()
    tmp = sorted(cols.keys())
    cols_rgb = [cvr.to_rgb(cols[k]) for k in tmp]
    intervals = np.array(tmp + [tmp[-1]+1]) - 0.5
    cmap, norm = colors.from_levels_and_colors(intervals,cols_rgb)
    
    data = np.zeros(shape=(args.N,args.L))

    fasta_sequences = SeqIO.parse(open(args.reads),'fasta')
    i = 0
    letters = ['A','C','G','T']
    for fasta in fasta_sequences:
        seq = str(fasta.seq).upper()
        #if seq.count('N')>0: continue
        for j in range(0,args.L):
            letter = seq[j]
            if letter not in letters:
                r = np.random.randint(4)
                letter = letters[r]
            if letter=='A': data[i,j] = 1
            elif letter=='C': data[i,j] = 2
            elif letter=='G': data[i,j] = 3
            else: data[i,j] = 4
        i += 1
        if i==args.N: break

        
    #np.random.shuffle(data)       
    #plotting
    plt.pcolor(data,cmap = cmap, norm = norm)
    plt.xlabel("sequence position")
    plt.suptitle(args.title)
    #plt.tight_layout()
    plt.savefig(args.outfile,dpi=300)
    plt.clf()
#end

seq2heatmap()
