'''
The purpose of this script is to combine two fastas for an AlphaFold screen

Input: 
    1. Path to a folder with all the fasta files for the proteins to screen against
    2. Path to a folder where the combined fasta files are saved

Output: 
    1. Folder with all one file with all fastas
'''

import os
import argparse

def is_fasta_file(file_name):
    return file_name.lower().endswith(".fasta")

def read_fasta_file(file_path, encoding='utf-8'):
    sequences = []
    try:
        with open(file_path, "r", encoding=encoding) as file:
            content = file.read()
            if '>' in content:
                sequences.append(content)
            else:
                print(f"Skipping file '{file_path}' as it does not contain FASTA format.")
    except Exception as e:
        print(f"Error processing file '{file_path}': {e}")
    
    return sequences

def write_fasta_file(file_path, sequences, encoding='utf-8'):
    with open(file_path, "w", encoding=encoding) as file:
        for sequence in sequences:
            file.write(sequence + "\n")

def main(fasta_folder_path, output_fasta_path):
    combined_sequences = []

    for filename in os.listdir(fasta_folder_path):
        if is_fasta_file(filename):
            input_file = os.path.join(fasta_folder_path, filename)
            try:
                sequences = read_fasta_file(input_file)
                combined_sequences.extend(sequences)
            except Exception as e:
                print(f"Error processing file '{filename}': {e}")

    write_fasta_file(output_fasta_path, combined_sequences)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine all FASTA files in a folder into one file.")
    parser.add_argument("fasta_folder_path", help="Path to the folder with all the FASTA files to combine.")
    parser.add_argument("output_fasta_path", help="Path to the output FASTA file.")
    
    args = parser.parse_args()
    
    main(args.fasta_folder_path, args.output_fasta_path)
