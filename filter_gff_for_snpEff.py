from Bio import SeqIO

with open('contigs_filtered_cdhit.fasta', 'r') as handle:
    contig_list = []
    for contig in SeqIO.parse(handle, 'fasta'):
        contig_list.append(contig.id)

with open('proteins_filtered_cdhit.fasta','r') as handle:
    protein_list = []
    for prot in SeqIO.parse(handle, 'fasta'):
        print(prot.id)
        protein_list.append(prot.id)

with open('Trinity.fasta.transdecoder.gff3', 'r') as handle:
    store_list = []
    for line in handle.readlines():
        if line.split('\t')[0] in contig_list:
                for prot in protein_list:
                    if line.find(prot) > -1:
                        store_list.append(line)

with open('Trinity.fasta.transdecoder.gff3_filtered', 'w') as handle:
    for line in store_list:
        handle.write(line)
     
