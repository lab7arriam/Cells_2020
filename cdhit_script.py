import subprocess 
import os
import argparse
from Bio import SeqIO

def write(seqs:list, file:str):
    with open(file, 'w') as handle:
        SeqIO.write(seqs, handle, 'fasta')

def read_fasta(file:str)->list:
    with open(file, 'r') as handle:
        seq_list = []
        for entry in SeqIO.parse(handle, 'fasta'):
            seq_list.append(entry)
    return seq_list

def cdhit_correction(contigs:str, proteins:str, cd_output:str, contig_output:str, threshold:int)->list:
    """
    first cluster proteins in proteins.fa with CD-Hit with a provided identity threshold
    then use the resulting file to /filter_empty/ onto the contigs.fa
    finally, write the filtered contig list to the file with /write/ function
    """
    DEVNULL = open(os.devnull, 'w')
    
    def is_tool(name)->bool:
        """
        Check whether `name` is on PATH
        """

        from distutils.spawn import find_executable

        return find_executable(name) is not None
    
    def cd_hit(proteins:str, cd_output:str, threshold:int):
        """
        run CD-Hit for protein sequences if the tool is in $PATH and write the output to the specified file
        """
    
        if is_tool('cd-hit'):
            subprocess.run(f'cd-hit -i {proteins} -o {cd_output} -c {threshold} -n 5 -d 0 -M 0 -T 72', shell=True,
                          stdout = DEVNULL, stderr = DEVNULL)
            
    def filter_missing(contigs:str, filtered_proteins:str):
        """
        wipe out contigs whose respective proteins were filtered on CD-Hit clustering stage
        """
        
        contig_list, protein_list = read_fasta(contigs), read_fasta(filtered_proteins)
        contig_list = (contig for contig in contig_list if contig.id in map(lambda x: x.id[:-3], protein_list))
        
        return contig_list
    
    cd_hit(proteins, cd_output, threshold)
    write(filter_missing(contigs, cd_output), contig_output)


if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('-nu', help='path to nucleotide file', type=str, required=True)
   parser.add_argument('-pr', help='path to protein file', type=str, required=True)
   parser.add_argument('-nuo', help='path to nucleotide output file', type=str, required=True)
   parser.add_argument('-pro', help='path to protein output file', type=str, required=True)
   parser.add_argument('-t', help='CD-Hit identity cutoff', type=float, required=True)

   args = parser.parse_args()
   nu, pr, nuo, pro, t = args.nu, args.pr, args.nuo, args.pro, args.t
   
   cdhit_correction(nu, pr, pro, nuo, t)
