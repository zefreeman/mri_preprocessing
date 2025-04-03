#!/bin/bash

# Define the base path to the derivatives folder
base_folder="/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/derivatives"

# Loop through all sub-IDXX/ses-mri01/func folders inside the derivatives folder
for func_folder in "$base_folder"/sub-ID*/ses-mri02/func; do
    if [ -d "$func_folder" ]; then
        # Loop through all files matching the pattern *task*_desc-confounds_timeseries.tsv
        for file in "$func_folder"/*task-TRMrun*_desc-preproc_bold.nii.gz; do
            echo "Processing file: $file"
            
            # Check if the file exists
            if [ -f "$file" ]; then
                # Read the header to find the index of the framewise_displacement column
                gunzip $file

	    else
                echo "Error: File $file not found"
            fi
        done
    else
        echo "Func folder not found in $sub_folder"
    fi
done