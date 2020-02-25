#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict

'''
Parse the clusters and assign sequence to the respective cluster name
'''

with open('proteins_joint_cdhit.fasta.clstr', 'r') as handle:
    clusters = defaultdict(list)
    for line in handle.readlines():
        if line.startswith('>Cluster'):
            root = line
        else:
            clusters[root].append(line.strip())

contigs_to_remove = []
contigs_to_append = []

for key in clusters.keys():
    if len(clusters[key]) > 1:
        for term in clusters[key]:
            if term.find('*') > -1 and term.find('comp') > -1:
               trinity_terms = [x for x in clusters[key] if x.find('TRINITY') > -1]
               champion = max(trinity_terms, key=lambda x: float(x.split(' ')[-1].rstrip('%')))
               print(champion)
               contigs_to_append.append(champion.split(' ')[-3].lstrip('>').rstrip('...'))
               contigs_to_remove.append(term.split(' ')[-2].lstrip('>').rstrip('...'))

with open('substitute_to_append_names.txt', 'w') as handle:
    for item in contigs_to_append:
        handle.write('{}\n'.format(item))

with open('substitute_to_remove_names.txt', 'w') as handle:
    for item in contigs_to_remove:
        handle.write('{}\n'.format(item))


with open('proteins_joint_cdhit.fasta', 'r') as handle:
    novel_proteins = []
    for protein in SeqIO.parse(handle, 'fasta'):
        if protein.id not in contigs_to_remove:
            novel_proteins.append(protein)

with open('proteins_joint.fasta', 'r') as handle:
    for protein in SeqIO.parse(handle, 'fasta'):
        if protein.id in contigs_to_append:
            novel_proteins.append(protein)

with open('contigs_joint.fasta', 'r') as handle:
    novel_contigs, protein_names = [], list(map(lambda x: x.id[:-3], novel_proteins))
    for contig in SeqIO.parse(handle, 'fasta'):
        if contig.id in protein_names:
            novel_contigs.append(contig)

SeqIO.write(novel_proteins, 'proteins_joint_cdhit_substituted.fasta', 'fasta')
SeqIO.write(novel_contigs, 'contigs_joint_cdhit_substituted.fasta', 'fasta')
