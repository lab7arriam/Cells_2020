#! /bin/bash
mkdir rcorrected
for dir in $(ls trimmed_loose); do if [[ $dir == 1 || $dir == 5 ]]; then
mkdir rcorrected/$dir
R1=$(ls trimmed_loose/$dir | grep 'R1')
R2=$(ls trimmed_loose/$dir | grep 'R2') 
perl ~/rcorrector/run_rcorrector.pl -1 trimmed_loose/$dir/$R1 -2 trimmed_loose/$dir/$R2 -od rcorrected/$dir -t 72
fi; done
