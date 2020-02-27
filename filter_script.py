#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import argparse

def read_fasta(file:str)->list:
    with open(file, 'r') as handle:
        seq_list = []
        for entry in SeqIO.parse(handle, 'fasta'):
            seq_list.append(entry)
    return seq_list


def filter_empty(contig_list:list, protein_list:list)->list:
    """
    wipe out contigs without any coding sequence embedded or coding sequence lost after CD-Hit clustering
    """
    contig_list = [contig for contig in contig_list if contig.id in map(lambda x: x.id[:-3], protein_list)]
    return contig_list


def filter_assembly(contig_list:list, protein_list:list)->list:
    """
    remove contigs which fail to contain a plausible coding sequence
    """
    
     def major_filter(contig_list:list, protein_list:list)->dict:
        """
        a major filtering function replacing filter_short_contig, filter_short_orfs and filte_overhangs functions
        """
        contig_dict = defaultdict(list)
        for contig in contig_list:
            if len(contig.seq) > 300: ## replaces filter_short_contigs in the first place
                for protein in protein_list:
                    if protein.id[:-3] == contig.id
                        if len(protein.seq)*3 >= len(contig.seq)*0.3 and len(protein.seq) >= 100 
                        ## replaces filter_short_orfs and filter_overhangs
                        contig_dict[contig.id].append(protein)
                    
        return contig_dict
        
    def leave_longest(contig_dict:dict)->dict:
        """
        leave only the longest CDS if more than one ORFs are found
        """
        contig_dict = dict(zip(contig_dict.keys(),
                           map(lambda x: max(x, key=lambda y:len(y.seq)), contig_dict.values())))
        return contig_dict    
    
    contig_list = filter_empty(contig_list, protein_list)
    contig_dict = major_filter(contig_list, protein_list)
    contig_dict = leave_longest(contig_dict)    
    contigs = (contig for contig in contig_list if contig.id in contig_dict.keys())
    
    return contigs, contig_dict.values()

def write(seqs:list, file:str):
    with open(file, 'w') as handle:
        SeqIO.write(seqs, handle, 'fasta')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('-pr', help='protein fasta file', type=str, required=True)
    parser.add_argument('-nu', help ='nucleotide fasta file', type=str, required=True)
    parser.add_argument('-o', help='output directory', type=str, required=True)

    args = parser.parse_args()
    pr, nu, o = args.pr, args.nu, args.o

    contig_list, protein_list = filter_assembly(read_fasta(nu), read_fasta(pr))
    write(contig_list, f'{o}/contigs_filtered.fasta')
    write(protein_list, f'{o}/proteins_filtered.fasta')

