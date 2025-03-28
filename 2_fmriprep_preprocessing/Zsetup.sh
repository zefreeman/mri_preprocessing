#  Author: Z Freeman
#    Date: March, 2023
# Purpose: Set up parameters for the MRIQC job and submit to cluster

# set shell to exit when a command returns an error (e)
# or when a variable is undefined (u) and print all 
# commands in stderr (x)
set -eux

# set up directory variables
DATA_DIR="/data/project/STAR/testBIDS/data" #
SCRIPTS_DIR="/data/project/STAR/testBIDS/scripts"
Q_LOG="${SCRIPTS_DIR}/q_logs/mriqc"
mkdir -p $Q_LOG 

# create input lists
sub_list="${SCRIPTS_DIR}/MRIQC_subject_list.txt"

# remove list if it already exists (Z: won't this remove the one we just made?)
if [[ -e $sub_list ]]; then rm $sub_list; fi

# find all subject directories in the bids folder
# and store them as an array
subj_in=( $( find $DATA_DIR/bids -type d -name "sub-*" ) )

# write each folder name in the input list
for i in ${subj_in[@]}; do echo $( basename $i ) >> $sub_list; done

# submit script
BOLD_ID=$( sbatch --parsable --array=1-$( wc -l < $sub_list )%25 \ 
-o $Q_LOG/bold_%A-%a.out $SCRIPTS_DIR/run_MRIQC_bold.job \  
$sub_list $DATA_DIR $SCRIPTS_DIR )

# check to see if the job is still running, and if it is
# then wait another 15 minutes before checking
while squeue | egrep "($BOLD_ID)"; do
    # job is in the queue
    echo "Still running MRIQC, waiting for 15 minutes..."
    sleep 15m
done
