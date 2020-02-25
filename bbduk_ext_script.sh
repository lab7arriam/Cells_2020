#! /bin/bash

abs_path="/home/yura/Ps_transcriptome/ext_2017"
for dir in $(ls $abs_path); do if [[ -d  $abs_path/$dir  ]]; then
/opt/bbmap/bbduk.sh ref=/opt/bbmap/resources/adapters.fa interleaved=f maq=20 ktrim=r  k=23  qtrim=rl t=72  minlength=50 in=$abs_path/$dir/$dir\_1.fastq  in2=$abs_path/$dir/$dir\_2.fastq  out=/home/yura/Ps_transcriptome/trimmed_ext/$dir/$dir\_R1_trimmed.fastq out2=/home/yura/Ps_transcriptome/trimmed_ext/$dir/$dir\_R2_trimmed.fastq overwrite=t
fi
done
