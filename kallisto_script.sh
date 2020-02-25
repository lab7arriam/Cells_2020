#! /bin/bash
for file in $(ls yeast_for_kallisto_3_rh | rev | cut -c 6- | rev| uniq )
do if
	short=$(echo $file | rev | cut -c 27- | rev)
	kallisto_linux-v0.44.0/kallisto quant -i indices/verified_orfs.idx -b 100 -o yeast_analized_5_rh/$short/ yeast_for_kallisto_2/$file\R1.fq yeast_for_kallisto_2/$file\R2.fq
   done
