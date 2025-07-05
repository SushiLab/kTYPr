#!/bin/bash

## Usage 
# fetch_paths.sh ./genomes/ ".*\.fna$" fetch.txt

# Inputs
search_dir=$1  # Directory to search in
regex=$2       # Regular expression to match
output_file=$3 # Output file to write matching files

# Ensure the search directory exists
if [[ ! -d "$search_dir" ]]; then
    echo "Directory does not exist: $search_dir"
    exit 1
fi

# Find files matching the regex and write to output file
find "$search_dir" -type f -regex "$regex" > "$output_file"

echo "Matching files written to: $output_file"