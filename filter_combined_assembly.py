#!/usr/bin/env python3

from Bio import SeqIO
import argparse


with open('proteins_filtered_cdhit.fasta', 'r') as handle:
    name_list = []
    for record in SeqIO.parse(handle, 'fasta'):
        name_list.append(record.id[:-3])

print(len(name_list))

with open('contigs_filtered.fasta', 'r') as handle:
    contig_list = []
    for record in SeqIO.parse(handle, 'fasta'):
        if record.id in name_list:
            contig_list.append(record)

SeqIO.write(contig_list, 'contigs_filtered_cdhit.fasta', 'fasta')

