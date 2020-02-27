#!/usr/bin/env python31

import os
from collections import defaultdict
import pandas as pd
import argparse

#wd = str(input())
#os.chdir(wd)

def exon_dict(e:str)->dict:
    """
    Stage 1: create an effective exon interval storage
    """
    # a SuperTranscript exon layout with extra column denoting a transcript of origin (preprocessed independently)
    gtf_default = pd.read_csv(e, sep='\t', header=None)
    gene_lib = defaultdict(list)
    exon_lib = defaultdict(list)
    for index, row in gtf_default.iterrows():
        row = list(row)
        if row[-1] not in gene_lib[row[0]]:
            gene_lib[row[0]].append(row[-1])
        exon_lib[row[-1]].append([row[3], row[4]])

    for key in exon_lib.keys():
        exon_lib[key].append(sum(map(lambda x: x[1] - x[0] + 1, exon_lib[key])))

    return gene_lib, exon_lib
    # it works!

def parse_VCF(VC:str)->list:
    """
    Stage 2: parse VCF, assess SNP location and identify transcript to bear the site in the final VCF
    If an exon refers to more than one transcript, create a new entry for each transcript
    """
    ### Parse VCF
    VC_file = str(input())
    vcf_default = []
    with open(VC_file, 'r') as handle:
        for i in handle.readlines():
            vcf_default.append(i.split('\t'))
    
    return vcf_default


### Create a list to store novel VC entries
def snp_trans_allocation(vcf_default:list, o:str):
    vcf_processed = []
    global gene_lib
    global exon_lib

    #### Determine supertranscript of origin
    for entry in vcf_default:
        for transcript in gene_lib[entry[0]]:
            entry_new = entry.copy()
            for i, exon in enumerate(exon_lib[transcript][:-1]):
                if exon[0] <= int(entry[1]) <= exon[1]:
                    entry_new[0] = transcript
                    ### Calculate the new position for the SNP site
                    new_pos = sum(map(lambda x: x[1] - x[0] + 1, exon_lib[transcript][:i])) + \
                          int(entry[1]) - exon_lib[transcript][i][0] + 1
                    entry_new[1] = str(new_pos)
                    vcf_processed.append('\t'.join(entry_new))
 
    output = o
    with open(output, 'w') as handle:
        for entry in vcf_processed:
            handle.write(entry)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', help='gtf file with exon layout', type=str, required=True)
    parser.add_argument('-vc', help='headless VCF file', type=str, required=True)
    parser.add_argument('-o', help='path to output file', type=str, required=True)

    args = parser.parse_args()
    e, vc, o = args.e, args.vc, args.o

    gene_lib, exon_lib = exon_dict(e)
    snp_trans_allocation(parse_VCF(vc), o)
