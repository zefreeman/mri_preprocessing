#!/bin/bash

set -eux

# set up directory variables
DATA_DIR="/data/project/STAR/second_test_fmriprep/data"
SCRIPTS_DIR="/data/project/STAR/second_test_fmriprep/scripts"
Q_LOG="${SCRIPTS_DIR}/q_logs/fmriprep"
#mkdir -p $Q_LOG

# create input lists
sub_list="${SCRIPTS_DIR}/fmriprep_subject_list.txt"

# remove list if it already exists
if [[ -e $sub_list ]]; then rm $sub_list; fi

# find all subject directories in the bids folder
# and store them as an array
subj_in=( $( find $DATA_DIR -type d -name "sub-*" ) )
echo "Found ${#subj_in[@]} directories matching pattern"

# write each folder name in the input list
for i in "${subj_in[@]}"; do 
  echo "$(basename "$i")" >> "$sub_list"
done
echo "Processed ${#subj_in[@]} directories"

# debug of sub_list unbound error
echo "sub_list is set to: ${sub_list}"

# submit script 
TESTPREP=$(sbatch --parsable --array=1-$( wc -l < $sub_list )%25 \
-o $Q_LOG/bold_%A-%a.out $SCRIPTS_DIR/Zfmriprep.job \
$sub_list $DATA_DIR $SCRIPTS_DIR )


#check to see if the job is still running, and if it is
#then wait another 15 minutes before checking
#while squeue | egrep "($TESTPREP)"; do
#    job is in the queue
#    echo "Still running fmriprep, waiting for 15 minutes..."
#    sleep 15m
#done
