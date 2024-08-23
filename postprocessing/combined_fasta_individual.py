'''
The purpose of this script is to combine two fastas for an AlphaFold screen

Input: 
    1. Path to a folder with all the fasta files for the proteins to screen against
    2. Path to fasta with the fasta file of your protein of interest
    3. Path to a folder where the combined fasta files are saved

Output: 
    1. Folder with all combined fasta files
'''

import os
import argparse

def is_fasta_file(file_name):
    return file_name.lower().endswith(".fasta")

def read_fasta_file(file_path, encoding='utf-8'):
    sequences = {}
    current_header = ""
    try:
        with open(file_path, "r", encoding=encoding) as file:
            content = file.read()

            # Check if the content contains the '>' symbol
            if '>' in content:
                lines = content.split('\n')
                for line in lines:
                    line = line.strip()
                    if line.startswith(">"):
                        current_header = line
                        sequences[current_header] = ""
                    else:
                        sequences[current_header] += line
            else:
                print(f"Skipping file '{file_path}' as it does not contain FASTA format.")
    except Exception as e:
        print(f"Error processing file '{file_path}': {e}")
        print(f"Contents of the file: {content}")

    return sequences

def write_fasta_file(file_path, sequences, encoding='utf-8'):
    with open(file_path, "w", encoding=encoding) as file:
        for header, sequence in sequences.items():
            file.write(header + "\n")
            file.write(sequence + "\n")

def main(fasta_folder_path, appended_fasta_path, output_fasta_path):
    for filename in os.listdir(fasta_folder_path):
        if is_fasta_file(filename):
            input_file = os.path.join(fasta_folder_path, filename)
            appended_fasta = appended_fasta_path

            try:
                existing_sequence = read_fasta_file(input_file)
                new_sequence = read_fasta_file(appended_fasta)
                existing_sequence.update(new_sequence)
                output_file = os.path.join(output_fasta_path, filename.replace(".fasta", "_PAX3_FOX01.fasta"))
                write_fasta_file(output_file, existing_sequence)
            except Exception as e:
                print(f"Error processing file '{filename}': {e}")
                print(f"Contents of the file: {read_file_content(input_file)}")

def read_file_content(file_path):
    with open(file_path, "r", encoding='utf-8') as file:
        return file.read()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine two FASTA files for AlphaFold screening.")
    parser.add_argument("fasta_folder_path", help="Path to the folder with all the FASTA files for the proteins to screen against.")
    parser.add_argument("appended_fasta_path", help="Path to the FASTA file of your protein of interest.")
    parser.add_argument("output_fasta_path", help="Path to the folder where the combined FASTA files will be saved.")
    
    args = parser.parse_args()
    
    main(args.fasta_folder_path, args.appended_fasta_path, args.output_fasta_path)
