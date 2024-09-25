# UniProt FASTA Fetcher

This repository contains a Python script and a Bash script to automate the process of fetching FASTA sequences from UniProt for a list of proteomes provided in a CSV file.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Python Script Usage](#python-script-usage)
  - [Bash Script Usage](#bash-script-usage)
- [Error Handling](#error-handling)
- [License](#license)
- [Contact](#contact)

## Requirements
- Python 3.x
- `requests` module (`pip install requests`)
- Bash (for batch processing with the provided script)

## Installation
1. Clone the repository:
   ```
   git clone https://github.com/Haayhur/MD_codes.git
   cd uniprot-fasta-fetcher
   ```
2. Install the required Python package:
   ```
   pip install requests
   ```

## Usage

### Python Script Usage

The Python script fetches a FASTA file for a specific UniProt Proteome ID from UniProt's REST API.

#### Command
```
python uniprot_fasta_fetcher.py <ProteomeID>
```
**ProteomeID**: The UniProt Proteome ID (e.g., `UP000005640`).

#### Example
```
python uniprot_fasta_fetcher.py UP000005640
```
This will download the FASTA file and save it as `UP000005640.fasta` in the current directory.

### Bash Script Usage

The Bash script automates batch processing of proteome IDs from a CSV file, fetching their corresponding FASTA files using the Python script.

#### Command
```
./fetch_fasta_from_csv.sh <csv_file>
```
**csv_file**: The CSV file containing proteome data. The CSV should be tab-separated and include a header. The first column should contain the proteome IDs.

#### Example CSV Structure
```
Proteome ID    Organism    Organism ID    Protein Count    BUSCO    CPD
UP000005640    Homo sapiens    9606        20350           95%      123
UP000006548    Escherichia coli 511145     4300            99%      456
```

#### Example Command
```
./fetch_fasta_from_csv.sh proteomes.csv
```
This will fetch the FASTA sequences for all proteome IDs in the CSV file and save them in the `fasta_files` directory.

## Error Handling
- **Python Script**: If the Python script fails to fetch a FASTA file, it will return an error message indicating the HTTP status code.
- **Bash Script**: The Bash script will output error messages if it cannot find the CSV file, Python script, or if a proteome ID fails to process.

## License
This project is licensed under the MIT License.

## Contact
For any questions or issues, please open an issue in the repository.
