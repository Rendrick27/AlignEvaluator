import os
import re


def read_fasta_files(folder_path):
    """
    Reads all FASTA files from the specified folder and
    returns a dictionary with species name and sequence.

    Args:
        folder_path (str): Path to the folder containing FASTA files.

    Returns:
        dict: A dictionary containing species names as keys and
         lists of sequences as values.
    """
    all_sequences = {}
    file_list = sorted(os.listdir(folder_path))  # Sorted list of file names
    for filename in file_list:
        if filename.endswith('.fasta'):
            file_path = os.path.join(folder_path, filename)
            sequences = read_fasta(file_path)
            for species, sequence in sequences.items():
                if species in all_sequences:
                    all_sequences[species].append(sequence)
                else:
                    all_sequences[species] = [sequence]
    return all_sequences


def read_fasta(file_path):
    """
    Reads a single FASTA file and returns a dictionary
    with species name and sequence.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary containing species names
        as keys and sequences as values.
    """
    sequences = {}
    with open(file_path, 'r') as file:
        current_species = None
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if current_species:
                    sequences[current_species] = sequence
                current_species = re.match(r'>(\w+)_\w+', line).group(1)
                sequence = ''
            else:
                sequence += line.strip()
        if current_species:
            sequences[current_species] = sequence
    return sequences


def write_combined_fasta(combined_sequences, output_file):
    """
    Writes combined sequences into a new FASTA file.

    Args:
        combined_sequences (dict): A dictionary containing species names
        as keys and concatenated sequences as values.
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, 'w') as file:
        for species, sequences in combined_sequences.items():
            sequence = ''.join(sequences)
            file.write(f'>{species}\n{sequence}\n')


# Folder containing the FASTA files
folder_path = 'dataset'

# Combine sequences from all FASTA files in the folder
combined_sequences = read_fasta_files(folder_path)

# Write the combined sequences into a new FASTA file
write_combined_fasta(combined_sequences, 'combined.fasta')

print("Combined sequences have been saved in 'combined.fasta'.")
