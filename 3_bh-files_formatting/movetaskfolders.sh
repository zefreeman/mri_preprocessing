#!/bin/bash

# Define the source and destination directories
source_dir="/Users/z/Documents/STAR_baseline_MRI/data/scans_N/"
destination_dir="/Users/z/Documents/STAR_baseline_MRI/data/scans_N/taskfiles/"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"
echo "Destination directory created: $destination_dir"

# Loop through all folders inside scans_L
for folder in "$source_dir"/*; do
    echo "Processing folder: $folder"
    if [ -d "$folder" ]; then
        echo "Folder is a directory: $folder"
        # Get the name of the folder
        folder_name=$(basename "$folder")
        echo "Folder name: $folder_name"

        # Check if the task folder exists within the subfolder "T1"
        task_folder="$folder/T1/task"
        echo "Task folder path: $task_folder"
        if [ -d "$task_folder" ]; then
            echo "Task folder exists: $task_folder"
            echo "Copying task folder from $task_folder to $destination_dir$folder_name"
            cp -r "$task_folder" "$destination_dir$folder_name"
            echo "Task folder copied successfully."
        else
            echo "Task folder does not exist in $folder"
        fi
    else
        echo "Skipping non-directory: $folder"
    fi
done

echo "Script execution completed."
