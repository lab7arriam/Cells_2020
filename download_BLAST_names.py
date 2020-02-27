#! /usr/bin/env python3

from Bio import Entrez as E
from Bio import SeqIO
import argparse
import sys

def parse_args(arg_list:list):
     parser = argparse.ArgumentParser()
     parser.add_argument('-i', type=str, required=True)
     parser.add_argument('-e', type=str, required=False, default='271296251017a@gmail.com')
     parser.add_argument('-db', type=str, required=False, default='protein')

     return parser.parse_args()
#     file, db = args.i, args.db





def read_file(file:str)->list:
    with open(file, 'r') as handle:
        return [x.strip() for x in handle.readlines()]


def fetch_names(names:list, db:str='protein')->list:
    
    out_list = []

    for i in range(0, len(names), 1000):
        query = ','.join(names[i:i+1000])
        extracted_names = SeqIO.parse(E.efetch(db=db, id=query, rettype='gb', retmode='txt'), 'gb')
        
        for j in range(i, min(i+1000, len(names))):
            if names[j] != 'None':
                out_list.append(next(extracted_names).description)
            else:
                out_list.append('NA')

    return out_list

def write_outputs(records:list, output:str):
    with open(output, 'w') as handle:
        for item in records:
            handle.write('%s\n' % item)


def main(args:list):
    pass

if __name__=='__main__':
     parser = argparse.ArgumentParser()
     parser.add_argument('-i', type=str, required=True)
     parser.add_argument('-db', type=str, required=False, default='protein')
     parser.add_argument('-o', type=str, required=True)

     args = parser.parse_args()
     file, db, output = args.i, args.db, args.o

     E.email, E.tool = 'yu.malovichko@arriam.ru', 'MyCustomScript'

     write_outputs(fetch_names(read_file(file), db=db), output)
     
