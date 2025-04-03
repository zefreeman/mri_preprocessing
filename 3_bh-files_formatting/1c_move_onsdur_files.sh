#!/bin/bash

# Source and destination folders
source_folder="/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/L2_onsdur_20250401/"
destination_folder="/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/derivatives/"

# Loop through all .mat files in the source folder
for file in "${source_folder}"*.mat; do
    # Extract sub-ID from the file name
    file_name=$(basename "$file")
    sub_id=$(echo "$file_name" | grep -oP '^sub-ID\d+')

    if [[ -n "$sub_id" ]]; then
        # Construct the destination directory
        destination_dir="${destination_folder}${sub_id}/ses-mri02/func/"

        # Ensure the destination directory exists, create it if necessary
        mkdir -p "$destination_dir"

        # Move the file to the destination directory
        mv "$file" "${destination_dir}${file_name}"
        echo "Moved $file_name to $destination_dir"
    else
        echo "Failed to find sub-ID in file: $file_name"
    fi
done