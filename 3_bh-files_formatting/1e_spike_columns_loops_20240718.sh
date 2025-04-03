#!/bin/bash

# Define the base path to the derivatives folder
base_folder="/data/project/STAR/STAR_preprocessed/L2_fmriprep/data/derivatives"

# Loop through all sub-IDXX/ses-mri01/func folders inside the derivatives folder
for func_folder in "$base_folder"/sub-ID*/ses-mri02/func; do
    if [ -d "$func_folder" ]; then
        # Loop through all files matching the pattern *task*_desc-confounds_timeseries.tsv
        for file in "$func_folder"/*task*_desc-confounds_timeseries.tsv; do
            echo "Processing file: $file"
            
            # Check if the file exists
            if [ -f "$file" ]; then
                # Read the header to find the index of the framewise_displacement column
                header=$(head -n 1 "$file")
                fw_index=$(echo "$header" | tr '\t' '\n' | grep -n 'framewise_displacement' | cut -d: -f1)
                
                # Check if framewise_displacement column exists in the file
                if [ -n "$fw_index" ]; then
                    # Find rows where framewise_displacement >= 0.5
                    rows=$(awk -v fw_index="$fw_index" -F'\t' 'NR>1 && $fw_index >= 0.5 {print NR}' "$file")
                    
                    # Process each row with framewise_displacement >= 0.5
                    while read -r row; do
                        # Create a new spike column for each row with zeros initially
                        col_name="spike_$((row-1))"
                        awk -v col_name="$col_name" 'BEGIN{FS=OFS="\t"} {print $0, (NR==1 ? col_name : 0)}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
                        
                        # Replace the zero with one in the corresponding row
                        awk -v row="$row" -F'\t' -v OFS='\t' 'NR==row {$(NF)=1} {print}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
                    done <<< "$rows"
                    
                    # Remove the spike_1 column if it exists
                    awk 'BEGIN{FS=OFS="\t"} {
                        if (NR==1) {
                            for (i=1; i<=NF; i++) {
                                if ($i == "spike_1") {
                                    col_to_remove = i;
                                }
                            }
                        }
                        if (col_to_remove) {
                            for (i=col_to_remove; i<NF; i++) {
                                $i = $(i+1);
                            }
                            NF--;
                        }
                        print
                    }' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
                    
                    echo "Processed file: $file"
                else
                    echo "Error: framewise_displacement column not found in $file"
                fi
            else
                echo "Error: File $file not found"
            fi
        done
    else
        echo "Func folder not found in $sub_folder"
    fi
done
