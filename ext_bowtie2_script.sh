#! /bin/bash

abs="/home/yura/Ps_transcriptome"
chp="trimmed_ext"
for dir in $(ls $abs/$chp); do if [[ -d $abs/$chp/$dir ]]; then
bowtie2 --threads 72 -x $abs/filtered_joined,BUSCOed/transcriptome_index -1 $abs/$chp/$dir/$dir\_R1_trimmed.fastq -2 $abs/$chp/$dir/$dir\_R2_trimmed.fastq -S $abs/ext_align_stats/$dir
fi
done
