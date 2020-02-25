#/bin/bash

abs='/home/yura/Ps_transcriptome/novel_merged_assembly'
for dir in $(ls $abs | grep 'filtered'); do
python3 /home/yura/BUSCO/busco-master/scripts/run_BUSCO.py -i `find $abs/$dir/* -name 'contigs*fasta'` -o $dir\.busco -l /home/yura/BUSCO/lineages/viridiplantae_odb10 -m tran &
done
