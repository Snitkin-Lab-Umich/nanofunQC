#!/bin/bash

# Check if the user has provided an argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <path_to_fastq_files>"
  exit 1
fi

# Get the path from the user input and remove any trailing slash
input_path=${1%/}

# Check if the provided path is a directory
if [ ! -d "$input_path" ]; then
  echo "Error: $input_path is not a valid directory."
  exit 1
fi

# Loop through all fastq.gz files in the provided directory
for file in "$input_path"/*.fastq.gz; do
  # Check if there are no fastq.gz files
  if [ ! -e "$file" ]; then
    echo "No fastq.gz files found in the directory."
    exit 1
  fi

  # Extract the sample name from the file name
  sample_name=$(basename "$file" .fastq.gz)
  
  # Create a directory with the sample name
  mkdir -p "$input_path/$sample_name"
  
  # Move the fastq.gz file into the created folder
  mv "$file" "$input_path/$sample_name/"
  
  # Check if the move was successful
  if [ $? -ne 0 ]; then
    echo "Error: Failed to move $file to $input_path/$sample_name"
    exit 1
  fi

  echo "Moved $file to $input_path/$sample_name"
done

echo "All fastq.gz files have been moved to their corresponding folders."
