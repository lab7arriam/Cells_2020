from Bio import SeqIO
import argparse

def read_fasta(file:str)->list:
    with open(file, 'r') as handle:
        seq_list = []
        for entry in SeqIO.parse(handle, 'fasta'):
            seq_list.append(entry)
    return seq_list

def find_by_length(seqs:list, cutoff:int):
    out = (seq for seq in seqs if len(seq.seq) >= cutoff)
    return out

def write(seqs:list, file:str):
    with open(file, 'w') as handle:
        SeqIO.write(seqs, handle, 'fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True)
    parser.add_argument('-o', type=str, required=True)
    parser.add_argument('-c', help='length cutoff', type=int, required=True)

    args = parser.parse_args()
    i, o, c = args.i, args.o, args.c
    
    write(find_by_length(read_fasta(i), c), o)

