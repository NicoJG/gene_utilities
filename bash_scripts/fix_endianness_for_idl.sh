#!/bin/bash

# Script to change ENDIANNESS = LITTLE to ENDIANNESS = BIG in output parameter files
# This is a specific fix for using the IDL on the IPP-HGW clusters
# Usage: ./fix_endianness_for_idl.sh file1 file2 ...
# Usage: ./fix_endianness_for_idl.sh *.txt
# Usage: ./fix_endianness_for_idl.sh "pattern*"
# Usage: ./fix_endianness_for_idl.sh /path/to/directory
# 
# If a directory is provided, it will recursively search for all files
# starting with "parameters" and apply the endianness change to them.

# Check if at least one argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename1|directory> [filename2] [...]"
    echo "Examples:"
    echo "  $0 file.txt"
    echo "  $0 *.txt"
    echo "  $0 file1.txt file2.txt"
    echo "  $0 /path/to/directory"
    echo ""
    echo "If a directory is provided, it will recursively search for all files"
    echo "starting with 'parameters' and apply the endianness change to them."
    exit 1
fi

# Initialize counters
files_processed=0
files_changed=0
total_replacements=0

# Function to process a single file
process_file() {
    local file="$1"
    local file_dir=$(dirname "$file")
    local idl_output_dir="$file_dir/idl_output"
    
    files_processed=$((files_processed + 1))
    
    # Check if the file contains any ENDIANNESS setting
    if grep -q "ENDIANNESS = \(LITTLE\|BIG\)" "$file"; then
        # Create idl_output directory if it doesn't exist
        if [ ! -d "$idl_output_dir" ]; then
            mkdir -p "$idl_output_dir"
            echo "Created directory: $idl_output_dir"
        fi
    fi
    
    # Check if the file contains the target string for replacement
    if grep -q "ENDIANNESS = LITTLE" "$file"; then
        # Count occurrences before replacement
        count=$(grep -c "ENDIANNESS = LITTLE" "$file")
        
        # Perform the replacement
        sed -i 's/ENDIANNESS = LITTLE/ENDIANNESS = BIG/g' "$file"
        
        echo "Changed: $file - Replaced $count occurrence(s) of 'ENDIANNESS = LITTLE' with 'ENDIANNESS = BIG'"
        files_changed=$((files_changed + 1))
        total_replacements=$((total_replacements + count))
    fi
}

# Process each argument
for arg in "$@"; do
    if [ -d "$arg" ]; then
        # If it's a directory, find all files starting with "parameters"
        echo "Searching directory: $arg"
        echo "Looking for files starting with 'parameters'..."
        
        # Find all files starting with "parameters" recursively
        while IFS= read -r -d '' file; do
            process_file "$file"
        done < <(find "$arg" -type f -name "parameters*" -print0)
        
    elif [ -f "$arg" ]; then
        # If it's a regular file, process it directly
        process_file "$arg"
        
    else
        # Handle wildcards and patterns
        found_any=false
        for file in $arg; do
            if [ -f "$file" ]; then
                process_file "$file"
                found_any=true
            fi
        done
        
        if [ "$found_any" = false ]; then
            echo "Warning: '$arg' is not a file or directory, or no matching files found"
        fi
    fi
done

# Summary
echo ""
echo "Summary:"
echo "  Files processed: $files_processed"
echo "  Files changed: $files_changed"
echo "  Total replacements: $total_replacements"