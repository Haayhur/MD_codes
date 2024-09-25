UniProt FASTA Fetcher
This repository contains a Python script and a Bash script to automate the process of fetching FASTA sequences from UniProt for a list of proteomes provided in a CSV file.

Table of Contents
Requirements
Installation
Usage
Python Script Usage
Bash Script Usage
Error Handling
License
Contact
Requirements
Python 3.x
requests module (pip install requests)
Bash (for batch processing with the provided script)
Installation
Clone the repository:

bash
코드 복사
git clone https://github.com/yourusername/uniprot-fasta-fetcher.git
cd uniprot-fasta-fetcher
Install the required Python package:

bash
코드 복사
pip install requests
Usage
Python Script Usage
The Python script fetches a FASTA file for a specific UniProt Proteome ID from UniProt's REST API.

Command
bash
코드 복사
python uniprot_fasta_fetcher.py <ProteomeID>
ProteomeID: The UniProt Proteome ID (e.g., UP000005640).
Example
bash
코드 복사
python uniprot_fasta_fetcher.py UP000005640
This will download the FASTA file and save it as UP000005640.fasta in the current directory.

Bash Script Usage
The Bash script automates batch processing of proteome IDs from a CSV file, fetching their corresponding FASTA files using the Python script.

Command
bash
코드 복사
./fetch_fasta_from_csv.sh <csv_file>
csv_file: The CSV file containing proteome data. The CSV should be tab-separated and include a header. The first column should contain the proteome IDs.
Example CSV Structure
yaml
코드 복사
Proteome ID    Organism    Organism ID    Protein Count    BUSCO    CPD
UP000005640    Homo sapiens    9606        20350           95%      123
UP000006548    Escherichia coli 511145     4300            99%      456
Example Command
bash
코드 복사
./fetch_fasta_from_csv.sh proteomes.csv
This will fetch the FASTA sequences for all proteome IDs in the CSV file and save them in the fasta_files directory.

Bash Script Details
The provided Bash script (fetch_fasta_from_csv.sh) works as follows:

It checks whether the CSV file and Python script are available.
It creates a directory called fasta_files to store the fetched FASTA files.
It reads the CSV file, skipping the header, and processes each proteome ID using the Python script.
It moves each fetched FASTA file to the fasta_files directory.
It prints the progress and status for each proteome ID processed.
bash
코드 복사
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
Error Handling
Python Script: If the Python script fails to fetch a FASTA file, it will return an error message indicating the HTTP status code.
Bash Script: The Bash script will output error messages if it cannot find the CSV file, Python script, or if a proteome ID fails to process.
License
This project is licensed under the MIT License.

Contact
For any questions or issues, feel free to reach out to [Your Name or Email].
