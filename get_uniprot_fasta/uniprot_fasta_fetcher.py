import requests
import sys

def fetch_uniprot_fasta(proteome_id):
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "fasta",
        "query": f"(proteome:{proteome_id})"
    }
    
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        return response.text
    else:
        return f"Error: Unable to fetch data. Status code: {response.status_code}"

def save_fasta(content, filename):
    with open(filename, 'w') as file:
        file.write(content)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <ProteomeID>")
        sys.exit(1)
    
    proteome_id = sys.argv[1]
    fasta_content = fetch_uniprot_fasta(proteome_id)
    
    if fasta_content.startswith("Error"):
        print(fasta_content)
    else:
        filename = f"{proteome_id}.fasta"
        save_fasta(fasta_content, filename)
        print(f"FASTA file saved as {filename}")
