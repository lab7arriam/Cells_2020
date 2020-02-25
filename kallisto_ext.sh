#! /bin/bash
for dir in $(ls trimmed_ext); do if [[ -d trimmed_ext/$dir ]]; then
	kallisto quant -i filtered_joined,BUSCOed/kalindex_31.idx -b 150 -o kallisto_ext/$dir trimmed_ext/$dir/$dir\_R1_trimmed.fastq trimmed_ext/$dir/$dir\_R2_trimmed.fastq
fi
done
