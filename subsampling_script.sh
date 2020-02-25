#! /bin/bash

abs='trimmed_ext'
dest='ext_subsamples'
for dir in $(ls $abs); do
END=4
mkdir $dest/$dir
for ((i=1;i<=END;i++)); do
    rand=$RANDOM
for file in $(ls  $abs/$dir); do
quant=$(awk 'END{printf("%.0f\n" , NR / 4 * 0.75)}' $abs/$dir/$file)
seqtk sample -s$rand $abs/$dir/$file $quant > $dest/$dir/$file\_sub$i 
done
done
done
