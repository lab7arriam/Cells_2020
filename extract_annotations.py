#!/usr/bin/env python3

import click
import os
import pandas as pd
from Bio import SeqIO

def store_annotations()->dict:
    annot_dict = dict()
    for i in os.listdir('.'):
        if i.find('.fa') > -1 and i.find('.p') == -1:
            with open(i, 'r') as handle:
                for entry in SeqIO.parse(handle, 'fasta'):
                    annot_dict[entry.id] = entry.description

    return annot_dict

@click.command()
@click.option('-t', help='A path to the BLAST outfmt6 table')
@click.option('-o', help='A path to the output table')

def parse_blast_output(t:str, o:str):
    annot_dict = store_annotations()
    input = pd.read_csv(t, sep='\t', header=None)
    output = pd.DataFrame(data={'Contig':['a']*len(input), 'Accession':['b']*len(input), 'Description':['c']*len(input)})
    for i in range(0, len(input)):
        key, contig = input.iloc[i,0], input.iloc[i,1]
        print(key, contig)
        output.iloc[i] = [contig, key, annot_dict[key]]

    output.to_csv(o, sep=';', header=True, index=False)

if __name__ == '__main__':
    parse_blast_output()
