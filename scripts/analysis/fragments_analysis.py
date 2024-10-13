from Bio import SeqIO
import os
import Levenshtein  # You need to install this package
from collections import Counter


def compare_fasta_sequences_with_edit_distance_and_nucleotide_count(fasta_files):
    # Dictionary to hold sequences from all files by sequence ID
    sequences_dict = {}

    # Loop through all FASTA files
    for fasta_file in fasta_files:
        # Parse each FASTA file
        counter = 0
        with open(fasta_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                counter+=1
                if counter == 50:
                    break
                seq_id = record.id
                sequence = str(record.seq)

                # If the sequence ID exists, append it to the dictionary entry
                if seq_id in sequences_dict:
                    sequences_dict[seq_id].append((fasta_file, sequence))
                else:
                    sequences_dict[seq_id] = [(fasta_file, sequence)]

    # Compare sequences with the same ID and calculate edit distances and nucleotide averages
    differences = []
    nucleotide_averages = {}
    print("finish  reading")
    for seq_id, file_seq_pairs in sequences_dict.items():
        first_file, first_seq = file_seq_pairs[0]  # Get the first sequence
        # Initialize counters for nucleotide counts
        nucleotide_counts = Counter()
        total_sequences = len(file_seq_pairs)

        for fasta_file, sequence in file_seq_pairs:
            # Count nucleotides in each sequence
            nucleotide_counts.update(sequence)

        # Calculate average nucleotide frequencies
        avg_nucleotides = {nuc: nucleotide_counts[nuc] / total_sequences for nuc in 'ATCG'}
        nucleotide_averages[seq_id] = avg_nucleotides

        # Compare the sequences and calculate edit distances
        for fasta_file, sequence in file_seq_pairs[1:]:
            if sequence != first_seq:
                # Calculate the edit (Levenshtein) distance
                edit_distance = Levenshtein.distance(first_seq, sequence)
                differences.append((seq_id, first_file, fasta_file, edit_distance))

    return differences, nucleotide_averages


# Example usage
fasta_files = ['1.fas', '2.fas', '3.fas']  # Add your 10 fasta file paths here

differences, nucleotide_averages = compare_fasta_sequences_with_edit_distance_and_nucleotide_count(fasta_files)
print(len(differences))
# Report the differences with edit distances
if differences:
    print(f"Found differences in sequences with the same name:")
    for seq_id, first_file, second_file, edit_distance in differences:
        print(f"Sequence ID: {seq_id}")
        print(f" - Files: {first_file} vs {second_file} | Edit Distance: {edit_distance}")
else:
    print("No differences found between sequences with the same names across the FASTA files.")

# Report the average nucleotide counts
print("\nAverage nucleotide counts for sequences:")
for seq_id, avg_nucleotides in nucleotide_averages.items():
    print(f"Sequence ID: {seq_id}")
    print(
        f" - A: {avg_nucleotides['A']:.2f}, T: {avg_nucleotides['T']:.2f}, C: {avg_nucleotides['C']:.2f}, G: {avg_nucleotides['G']:.2f}")
