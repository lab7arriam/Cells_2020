#!/usr/bin/env python3

from Bio import SeqIO

with open('missing_proteins.txt', 'r') as handle:
    miss_list = []
    for line in handle.readlines():
        miss_list.append(line.rstrip('\n'))

with open('contigs_filtered.fasta', 'r') as handle:
    contig_list = []
    for contig in SeqIO.parse(handle, 'fasta'):
        #print(contig.id)
        if contig.id in miss_list:
            contig_list.append(contig)

with open('proteins_filtered.fasta', 'r') as handle:
    protein_list = []
    for protein in SeqIO.parse(handle, 'fasta'):
        if protein.id[:-3] in miss_list:
            protein_list.append(protein)

SeqIO.write(contig_list, 'missing_contigs.fasta', 'fasta')
SeqIO.write(protein_list, 'missing_proteins.fasta', 'fasta')
