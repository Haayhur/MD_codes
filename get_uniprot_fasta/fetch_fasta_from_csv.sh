#!/bin/bash

# Check if the CSV file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

CSV_FILE=$1

# Check if the CSV file exists
if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file not found: $CSV_FILE"
    exit 1
fi

# Check if the Python script exists
PYTHON_SCRIPT="uniprot_fasta_fetcher.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create a directory for output files
OUTPUT_DIR="fasta_files"
mkdir -p "$OUTPUT_DIR"

# Read the CSV file, skipping the header
tail -n +2 "$CSV_FILE" | while IFS=$'\t' read -r proteome_id organism organism_id protein_count busco cpd
do
    echo "Processing Proteome ID: $proteome_id"
    
    # Run the Python script
    python "$PYTHON_SCRIPT" "$proteome_id"
    
    # Move the resulting FASTA file to the output directory
    mv "${proteome_id}.fasta" "$OUTPUT_DIR/" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        echo "FASTA file saved in $OUTPUT_DIR/${proteome_id}.fasta"
    else
        echo "Error: Failed to process $proteome_id"
    fi
    
    echo "------------------------"
done

echo "Batch processing complete. FASTA files are stored in $OUTPUT_DIR/"
