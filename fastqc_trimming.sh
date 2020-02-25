#! /bin/bash
for dir in $(ls trimmed_loose)
do
	sudo docker run -v /home/yura/Ps_transcriptome:/yu -w /yu --rm quay.io/biocontainers/fastqc:0.11.7--5 fastqc -t 72 trimmed_loose/$dir/$dir\_R1_trimmed.fastq trimmed_loose/$dir/$dir\_R2_trimmed.fastq -o=/yu/fastqc_reports_trimmed
done
