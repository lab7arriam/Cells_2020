#!/usr/bin/env python3

from Bio import SeqIO
import argparse

def read_fasta(file:str)->list:
    with open(file, 'r') as handle:
        seq_list = []
        for entry in SeqIO.parse(handle, 'fasta'):
            seq_list.append(entry)
    return seq_list

def mean_length(seqs:list)->str:
    """
    return the longest, the shortest and the mean lengths in the fasta file
    """
    shortest = len(min(seqs, key = lambda x: len(x.seq)).seq)
    longest = len(max(seqs, key = lambda x: len(x.seq)).seq)
    mean = sum(list(map(lambda x: len(x.seq), seqs))) / len(seqs)

    return f'The longest sequence is {longest} symbols long, the shortest sequence is {shortest} symbols long, mean sequence length across the file is {mean} symbols'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='input file path', type=str, required=True)
    parser.add_argument('-q', help='list of contigs to be processed', type=str, required=True)

    args = parser.parse_args()
    i, q = args.i, args.q
    with open(q, 'r') as handle:
        query_list = []
        for line in handle.readlines():
            query_list.append(line.rstrip('\n'))
    print(mean_length([x for x in read_fasta(i) if x.id in query_list]))
