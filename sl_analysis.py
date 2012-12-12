#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Spliced-leader analysis
Keith Hughitt <khughitt@umd.ed>
2012/11/26

Analyzes an input file containing information about RNA-seq contigs for which
spliced leader (SL) reads were mapped inside of coding regions of the reference
genome. 

The script searches for the next 1-2 start codons following the SL sequence,..
"""
from Bio import SeqIO
from pandas import read_csv, Series

# Input files
ref_genome = '../data/TcruziEsmeraldo-LikeGenomic_TriTrypDB-4.1.fasta'
problematic_sl_contigs = '../data/inside_CDS_total_5UTR_coodinates_sequences.scoring.txt'

# Read in T. cruzi reference genome
seqs = dict((seq.id[7:-2], seq) for seq in SeqIO.parse(ref_genome, 'fasta'))

# Parse ambiguous problematic trans-splicing table
cols = ['gene id', 'strand', 'chromosome',
        'trans_splicing_start', 'cds_end', 'frequency', 'score']

data = read_csv(problematic_sl_contigs, sep=' ', index_col=0, 
                usecols=[0, 1, 2, 3, 4, 6, 8], names=cols)

# Find start sites
start_sites = []

for i, gene in enumerate(data.values):
    # get gene sequence
    #[-, Chr24, 249358, 251338, 1, 1]
    chromosome = seqs["Tc" + gene[1]]
    bioseq = chromosome[gene[2]:gene[3]]

    # if on negative strand, take reverse complement
    if gene[0] == '-':
        bioseq = bioseq.reverse_complement()

    # return next start-site if one exists
    start_sites.append(bioseq.seq.find('ATG'))

# add start site column to right of existing data
data = data.join(Series(start_sites, index=data.index, name='ATG'))
data.to_csv('../data/sl_analysis_output.csv')


