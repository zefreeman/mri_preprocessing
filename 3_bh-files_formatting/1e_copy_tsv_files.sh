#!/bin/bash

# Define the base path to the derivatives folder
base_folder="/data/project/STAR/STAR_preprocessed/N_fmriprep/data/derivatives"

# Loop through all sub-IDXX folders inside the derivatives folder
for sub_folder in "$base_folder"/sub-ID*/ses-mri01; do
    if [ -d "$sub_folder" ]; then
        # Define the original_tsv folder path
        original_tsv_folder="$sub_folder/original_tsv"
        
        # Create the original_tsv folder if it doesn't exist
        mkdir -p "$original_tsv_folder"
        
        # Define the path to the func folder
        func_folder="$sub_folder/func"
        
        # Check if the func folder exists
        if [ -d "$func_folder" ]; then
            # Loop through all files matching the pattern *task*_desc-confounds_timeseries.tsv
            for file in "$func_folder"/*task*_desc-confounds_timeseries.tsv; do
                if [ -f "$file" ]; then
                    echo "Copying file: $file to $original_tsv_folder"
                    # Copy the file to the original_tsv folder
                    cp -i "$file" "$original_tsv_folder"
                else
                    echo "No matching files found in $func_folder"
                fi
            done
        else
            echo "Func folder not found in $sub_folder"
        fi
    else
        echo "Ses-mri01 folder not found in $sub_folder"
    fi
done
