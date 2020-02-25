#! /bin/bash

for dir in $(ls ~/Ps_transcriptome/trimmed_loose);do
if [[ $dir -ne 1 ]] && [[ $dir -ne 5 ]]; then
nohup kallisto quant -i for_kallisto.idx -b 200 -o ~/Ps_transcriptome/kallisto_results/$dir ../trimmed_loose/$dir/$dir\_R1_trimmed.fastq ../trimmed_loose/$dir/$dir\_R2_trimmed.fastq 2>$dir\_err.txt &
fi
done
