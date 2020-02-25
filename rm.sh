#! /bin/bash
abs='/home/yura/Ps_transcriptome'
tl='trimmed_loose'
for dir in $(ls $abs/$tl);
do
for file in $(ls $abs/$tl/$dir);
do
/opt/bbmap/removemicrobes.sh in=$abs/$tl/$dir/$file outu=$abs/$tl/$dir\_rm/$file\_rm.fq -outm=$abs/$tl/$dir\_trashrm/$file\_trash threads=72 build=3;
done;
done

