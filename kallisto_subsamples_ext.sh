#! /bin/bash
abs='/home/yura/Ps_transcriptome/ext_subsamples'
for dir in $(ls $abs); do if [[ -d $abs/$dir ]]; then
num=4
for ((i=1;i<=num;i++)); do
	kallisto quant -i trepnr_cdhit_09.idx -b 150 -o kallisto_ext_TE/$dir\_$i $abs/$dir/$dir\_R1_trimmed.fastq_sub$i $abs/$dir/$dir\_R2_trimmed.fastq_sub$i
done
fi
done
