<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>UniProt FASTA Fetcher</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
        }
        h1, h2, h3 {
            color: #333;
        }
        code {
            background-color: #f4f4f4;
            padding: 2px 5px;
            border-radius: 4px;
            font-size: 1rem;
        }
        pre {
            background-color: #f4f4f4;
            padding: 10px;
            border-radius: 5px;
            overflow-x: auto;
        }
        a {
            color: #007bff;
            text-decoration: none;
        }
    </style>
</head>
<body>

<h1>UniProt FASTA Fetcher</h1>

<p>This repository contains a Python script and a Bash script to automate the process of fetching FASTA sequences from UniProt for a list of proteomes provided in a CSV file.</p>

<h2>Table of Contents</h2>
<ul>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a>
        <ul>
            <li><a href="#python-script-usage">Python Script Usage</a></li>
            <li><a href="#bash-script-usage">Bash Script Usage</a></li>
        </ul>
    </li>
    <li><a href="#error-handling">Error Handling</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
</ul>

<h2 id="requirements">Requirements</h2>
<ul>
    <li>Python 3.x</li>
    <li><code>requests</code> module (<code>pip install requests</code>)</li>
    <li>Bash (for batch processing with the provided script)</li>
</ul>

<h2 id="installation">Installation</h2>
<ol>
    <li>Clone the repository:
        <pre><code>git clone https://github.com/yourusername/uniprot-fasta-fetcher.git
cd uniprot-fasta-fetcher</code></pre>
    </li>
    <li>Install the required Python package:
        <pre><code>pip install requests</code></pre>
    </li>
</ol>

<h2 id="usage">Usage</h2>

<h3 id="python-script-usage">Python Script Usage</h3>

<p>The Python script fetches a FASTA file for a specific UniProt Proteome ID from UniProt's REST API.</p>

<h4>Command</h4>
<pre><code>python uniprot_fasta_fetcher.py &lt;ProteomeID&gt;</code></pre>
<p><strong>ProteomeID</strong>: The UniProt Proteome ID (e.g., <code>UP000005640</code>).</p>

<h4>Example</h4>
<pre><code>python uniprot_fasta_fetcher.py UP000005640</code></pre>
<p>This will download the FASTA file and save it as <code>UP000005640.fasta</code> in the current directory.</p>

<h3 id="bash-script-usage">Bash Script Usage</h3>

<p>The Bash script automates batch processing of proteome IDs from a CSV file, fetching their corresponding FASTA files using the Python script.</p>

<h4>Command</h4>
<pre><code>./fetch_fasta_from_csv.sh &lt;csv_file&gt;</code></pre>
<p><strong>csv_file</strong>: The CSV file containing proteome data. The CSV should be tab-separated and include a header. The first column should contain the proteome IDs.</p>

<h4>Example CSV Structure</h4>
<pre><code>Proteome ID    Organism    Organism ID    Protein Count    BUSCO    CPD
UP000005640    Homo sapiens    9606        20350           95%      123
UP000006548    Escherichia coli 511145     4300            99%      456</code></pre>

<h4>Example Command</h4>
<pre><code>./fetch_fasta_from_csv.sh proteomes.csv</code></pre>
<p>This will fetch the FASTA sequences for all proteome IDs in the CSV file and save them in the <code>fasta_files</code> directory.</p>

<h3>Bash Script Details</h3>

<p>The provided Bash script (<code>fetch_fasta_from_csv.sh</code>) works as follows:</p>
<ol>
    <li>It checks whether the CSV file and Python script are available.</li>
    <li>It creates a directory called <code>fasta_files</code> to store the fetched FASTA files.</li>
    <li>It reads the CSV file, skipping the header, and processes each proteome ID using the Python script.</li>
    <li>It moves each fetched FASTA file to the <code>fasta_files</code> directory.</li>
    <li>It prints the progress and status for each proteome ID processed.</li>
</ol>

<pre><code>#!/bin/bash

# Check if the CSV file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 &lt;csv_file&gt;"
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
    mv "${proteome_id}.fasta" "$OUTPUT_DIR/" 2&gt;/dev/null
    
    if [ $? -eq 0 ]; then
        echo "FASTA file saved in $OUTPUT_DIR/${proteome_id}.fasta"
    else
        echo "Error: Failed to process $proteome_id"
    fi
    
    echo "------------------------"
done

echo "Batch processing complete. FASTA files are stored in $OUTPUT_DIR/"</code></pre>

<h2 id="error-handling">Error Handling</h2>
<ul>
    <li><strong>Python Script</strong>: If the Python script fails to fetch a FASTA file, it will return an error message indicating the HTTP status code.</li>
    <li><strong>Bash Script</strong>: The Bash script will output error messages if it cannot find the CSV file, Python script, or if a proteome ID fails to process.</li>
</ul>

<h2 id="license">License</h2>
<p>This project is licensed under the MIT License.</p>

<h2 id="contact">Contact</h2>
<p>For any questions or issues, feel free to reach out to <a href="mailto:youremail@example.com">Your Name</a>.</p>

</body>
</html>
