'''
The purpose of this script is to extract all the uniprot ids from the input fasta file and output a txt file with just the uniprot ids

Input:
    1. fasta_file_path - Path to the FASTA file to extract the Uniprot IDs from
    2. output_file_path - Path to the output text file to save the Uniprot IDs. **Make sure you specifiy the output file to be a txt file**
'''
import re
import argparse

def extract_uniprot_ids(fasta_file_path):
    uniprot_ids = []
    with open(fasta_file_path, "r", encoding='utf-8') as file:
        for line in file:
            if line.startswith(">"):
                # Extract the UniProt ID using regex
                match = re.search(r"\|(\w+)\|", line)
                if match:
                    uniprot_id = match.group(1)
                    uniprot_ids.append(uniprot_id)
    return uniprot_ids

def save_uniprot_ids(uniprot_ids, output_file_path):
    with open(output_file_path, "w", encoding='utf-8') as file:
        for uniprot_id in uniprot_ids:
            file.write(uniprot_id + "\n")

def main(fasta_file_path, output_file_path):
    uniprot_ids = extract_uniprot_ids(fasta_file_path)
    save_uniprot_ids(uniprot_ids, output_file_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract UniProt IDs from a FASTA file and save them to a text file.")
    parser.add_argument("fasta_file_path", help="Path to the FASTA file to extract UniProt IDs from.")
    parser.add_argument("output_file_path", help="Path to the output text file to save the UniProt IDs.")
    
    args = parser.parse_args()
    
    main(args.fasta_file_path, args.output_file_path)
