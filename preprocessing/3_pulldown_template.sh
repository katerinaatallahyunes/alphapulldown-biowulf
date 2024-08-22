#!/bin/sh
mkdir folder_name   
#Count the number of jobs corresponding to the number of sequences:
baits=`grep -c "" your_file.txt` #count lines even if the last one has no end of line
candidates=`grep -c "" your_file.txt` #count lines even if the last one has no end of line
count=$(( $baits * $candidates ))
#Run the job array, 10 jobs at a time:
sbatch --array=1-$count%10 2_pulldown_template.sh
