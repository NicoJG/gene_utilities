#!/bin/bash

# Script to combine all scan.log files from subdirectories into one big file
# Usage: ./summarize_scan_logs.sh [input_directory] [output_file]
# 
# Safety features:
# - Prompts before overwriting existing files
# - Validates input directory exists
# - Provides clear error messages

# Default values
INPUT_DIR="${1:-.}"  # Current directory if not specified
OUTPUT_FILE="${2:-combined_scan_logs.txt}"  # Default output filename

# Validate arguments
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    echo "Usage: $0 [input_directory] [output_file]"
    echo ""
    echo "Combines all scan.log files found in subdirectories into one file."
    echo ""
    echo "Arguments:"
    echo "  input_directory   Directory to search for scan.log files (default: current directory)"
    echo "  output_file       Output filename (default: combined_scan_logs.txt)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/gene/runs my_scan_summary.txt"
    exit 0
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' does not exist."
    echo "Please provide a valid directory path."
    echo "Usage: $0 [input_directory] [output_file]"
    exit 1
fi

# Convert to absolute path for clearer output
INPUT_DIR=$(realpath "$INPUT_DIR")

# Check if we have read permissions
if [ ! -r "$INPUT_DIR" ]; then
    echo "Error: No read permission for directory '$INPUT_DIR'."
    exit 1
fi

# Ensure output file has a safe name (not a directory)
if [ -d "$OUTPUT_FILE" ]; then
    echo "Error: '$OUTPUT_FILE' is a directory, not a file."
    echo "Please specify a filename for the output."
    exit 1
fi

# Check if output file already exists and ask for confirmation
if [ -f "$OUTPUT_FILE" ]; then
    echo "Warning: Output file '$OUTPUT_FILE' already exists."
    echo -n "Do you want to overwrite it? (y/N): "
    read -r response
    case "$response" in
        [yY]|[yY][eE][sS])
            echo "Overwriting existing file: $OUTPUT_FILE"
            rm "$OUTPUT_FILE"
            ;;
        *)
            echo "Operation cancelled. Please choose a different output filename."
            echo "Usage: $0 [input_directory] [output_file]"
            exit 1
            ;;
    esac
fi

echo "Searching for scan.log files in: $INPUT_DIR"
echo "Output file: $OUTPUT_FILE"
echo "Starting at: $(date)"
echo

# Initialize counter
count=0

# Find all scan.log files and process them
while IFS= read -r -d '' logfile; do
    # Get the directory path of the scan.log file
    logdir=$(dirname "$logfile")
    
    echo "Processing: $logfile"
    
    # Add header with path information
    {
        echo "================================================================================"
        echo "SCAN LOG FILE: $logfile"
        #echo "DIRECTORY: $logdir"
        #echo "TIMESTAMP: $(date)"
        echo "================================================================================"
        echo
    } >> "$OUTPUT_FILE"
    
    # Append the content of the scan.log file
    cat "$logfile" >> "$OUTPUT_FILE"
    
    # Add separator
    {
        echo
        echo
        echo "################################################################################"
        echo
    } >> "$OUTPUT_FILE"
    
    ((count++))
    
done < <(find "$INPUT_DIR" -name "scan.log" -type f -print0)

# Summary
echo
echo "Summary:"
echo "========="
echo "Found and processed $count scan.log files"
echo "Combined output saved to: $OUTPUT_FILE"
echo "Total lines in output: $(wc -l < "$OUTPUT_FILE" 2>/dev/null || echo "0")"
echo "Finished at: $(date)"

if [ $count -eq 0 ]; then
    echo
    echo "Warning: No scan.log files found in '$INPUT_DIR' or its subdirectories."
    exit 1
fi
