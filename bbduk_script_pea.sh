#! /bin/bash

abs_path="/home/yura/Ps_transcriptome/072017_pisum_transcriptome"
for dir in $(ls $abs_path); do if [[ -d  $abs_path/$dir  ]] && [[ $dir > 6 ]]; then for file in $(ls $abs_path/$dir); do
/opt/bbmap/bbduk.sh ref=/opt/bbmap/resources/adapters.fa interleaved=f maq=20 ktrim=r  k=23  qtrim=rl t=72  minlength=50 in=$abs_path/$dir/$dir\_1.fastq.gz  in2=$abs_path/$dir/$dir\_2.fastq.gz  out=/home/yura/Ps_transcriptome/trimmed_loose/$dir/$dir\_R1_trimmed.fastq out2=/home/yura/Ps_transcriptome/trimmed_loose/$dir/$dir\_R2_trimmed.fastq overwrite=t
done
fi
done
