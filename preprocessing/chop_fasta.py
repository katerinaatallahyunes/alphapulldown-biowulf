'''
The purpose of this script is to chop full length fastas into shorter sequences

Input: 
    1. Full length amino acid sequence fasta file (User-input path)
    2. Start and end residue numbers of amino acid sequence 
Output: Segment of full length amino acid sequence fasta file (Saved to output 
'''
import argparse

'''
Reads a FASTA file and returns a dictionary where keys are sequence IDs and values are sequences.
'''
def read_fasta(file_path):
    sequences = {}
    current_id = None
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            elif current_id is not None:
                sequences[current_id] += line
    return sequences

'''
Extracts subsequence from a given sequence based on start and end indices.
'''
def extract_subsequence(sequence, start_index, end_index):
    return sequence[start_index - 1:end_index]  # Adjust for 1-based indexing

def main(input_file, output_file, start_index, end_index):
    # Read the FASTA file
    fasta_sequences = read_fasta(input_file)

    # Check if there is at least one sequence in the file
    if not fasta_sequences:
        print("No sequences found in the FASTA file.")
        return

    # Extract the first sequence and its ID
    first_sequence_id, first_sequence = list(fasta_sequences.items())[0]
    first_sequence_id = first_sequence_id + " " + str(start_index) + "-" + str(end_index)
    print(first_sequence_id)

    # Extract the subsequence
    subsequence = extract_subsequence(first_sequence, start_index, end_index)

    # Write the subsequence to the output file
    with open(output_file, "w") as output:
        output.write(f">{first_sequence_id}\n")
        output.write(subsequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Chop full length fastas into shorter segments based on sequence.")
    parser.add_argument("input_file", help="Path to the input FASTA file.")
    parser.add_argument("output_file", help="Path to save the output segment of the sequence.")
    parser.add_argument("start_index", type=int, help="Start residue number of the amino acid sequence.")
    parser.add_argument("end_index", type=int, help="End residue number of the amino acid sequence.")

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.start_index, args.end_index)