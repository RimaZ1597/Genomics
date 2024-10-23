"""
Task : assemble small reads to contigs with greedy approch

"""


import pandas as pd
import matplotlib.pyplot as plt


def find_overlap(read1, read2):
    """
    Calculate the overlap between two reads.
    """
    max_overlap = 0
    overlap_sequence = ""
    min_len = min(len(read1), len(read2))

    for i in range(1, min_len + 1):
        if read1[-i:] == read2[:i]:
            max_overlap = i
            overlap_sequence = read1[-i:]
    return max_overlap, overlap_sequence


def merge_reads_with_overlap(read1, read2, overlap_len):
    """
    Merge two reads based on the overlap length.
    """
    return read1 + read2[overlap_len:]


def write_contigs_to_fasta(contigs, output_file):
    """
    Write contigs to a FASTA file.
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for i, contig in enumerate(contigs):
            f.write(f">Contig_{i + 1}\n{contig}\n")


def assemble_reads(reads_list, k):
    """
    Assemble reads into contigs based on overlaps.
    """
    contigs = reads_list[:]  # Make a copy of the list of reads
    final_contigs = []

    while contigs:
        current_read = contigs.pop(0)
        merged = True

        while merged:
            merged = False
            for i, next_read in enumerate(contigs):
                overlap_len, overlap_sequence = find_overlap(current_read, next_read)

                if overlap_len >= k:
                    merged_contig = merge_reads_with_overlap(current_read, next_read, overlap_len)

                    print(f"Merging:\nRead1: {current_read}\nRead2: {next_read}\n"
                          f"Overlap Length: {overlap_len}, Overlap Sequence: {overlap_sequence}\n"
                          f"Resulting Contig: {merged_contig}\n")

                    current_read = merged_contig
                    del contigs[i]
                    merged = True
                    break  # Restart with the updated contig

        final_contigs.append(current_read)

    return final_contigs


# Load the data
df = pd.read_csv("seqReadFile2023.txt", header=None, names=["Read"])
reads_list = df["Read"].tolist()

K_VALUE = 10  # Set the overlap length parameter

# Assemble the reads
final_assembled_contigs = assemble_reads(reads_list, K_VALUE)

# Save the final contigs to a CSV file
OUTPUT_CSV_FILE = "merged_all_contigs.csv"
contigs_df = pd.DataFrame(final_assembled_contigs, columns=["Contig"])
contigs_df.to_csv(OUTPUT_CSV_FILE, index=False)
print(f"Final merged contigs saved to '{OUTPUT_CSV_FILE}'")


def find_contig_overlaps(contigs):
    """
    Find overlaps between contigs.
    """
    overlaps = []
    for i in range(len(contigs)):
        for j in range(i + 1, len(contigs)):
            if contigs[i] in contigs[j] or contigs[j] in contigs[i]:
                overlaps.append((contigs[i], contigs[j]))
    return overlaps


def merge_overlapping_contigs(overlapping_contigs):
    """
    Merge overlapping contigs into a single representation.
    """
    merged = {}

    for contig1, contig2 in overlapping_contigs:
        if contig1 in merged:
            merged[contig1] = merged[contig1]
        elif contig2 in merged:
            merged[contig2] = merged[contig2]
        else:
            if contig1 in contig2:
                merged[contig2] = contig2
            elif contig2 in contig1:
                merged[contig1] = contig1
            else:
                merged[contig1] = contig1
                merged[contig2] = contig1 + contig2[len(contig1):]

    return list(merged.values())


# Find overlaps in the contigs
overlapping_contigs_list = find_contig_overlaps(final_assembled_contigs)

# Merge any overlaps found
merged_contigs_list = merge_overlapping_contigs(overlapping_contigs_list)

# Check for overlaps in the merged contigs
new_overlapping_contigs_list = find_contig_overlaps(merged_contigs_list)

# Output the results for new overlaps
if new_overlapping_contigs_list:
    print("Found overlaps in merged contigs:")
    for pair in new_overlapping_contigs_list:
        print(f"{pair[0]} \n<->\n {pair[1]}\n")

    merged_contigs_list = merge_overlapping_contigs(new_overlapping_contigs_list)

    # Print the final merged contigs
    print("Final merged contigs:")
    for contig in merged_contigs_list:
        print(contig)
else:
    print("No overlaps found in merged contigs.")

# Check for overlaps in the final merged contigs
final_overlapping_contigs_list = find_contig_overlaps(merged_contigs_list)

# Output the results for final overlaps
if final_overlapping_contigs_list:
    print("Found overlaps in final merged contigs:")
    for pair in final_overlapping_contigs_list:
        print(f"{pair[0]} \n<->\n {pair[1]}\n")
else:
    print("No overlaps found in final merged contigs.")

# Merge any overlaps found in final merged contigs
if final_overlapping_contigs_list:
    merged_final_contigs_list = merge_overlapping_contigs(final_overlapping_contigs_list)

    # Print the newly merged final contigs
    print("Merged final contigs after resolving overlaps:")
    for contig in merged_final_contigs_list:
        print(contig)

    # Write the final merged contigs to a FASTA file again if necessary
    FINAL_FASTA_OUTPUT = "RimaAssignment1.fasta"
    write_contigs_to_fasta(merged_final_contigs_list, FINAL_FASTA_OUTPUT)
    print(f"Final merged contigs saved to FASTA format in '{FINAL_FASTA_OUTPUT}'")
else:
    print("No final merged contigs to save.")


def calculate_coverage(all_reads, merged_contigs):
    """
    Calculate coverage of reads across merged contigs.
    """
    coverage = {contig: 0 for contig in merged_contigs}

    for read in all_reads:
        for contig in merged_contigs:
            if read in contig:
                coverage[contig] += 1

    return coverage


def plot_coverage(coverage_data):
    """
    Plot coverage of reads across merged contigs using generic labels (contig1, contig2, etc.).
    """
    contig_labels = [f"contig{i + 1}" for i in range(len(coverage_data))]  # Create generic contig labels
    coverage_counts = list(coverage_data.values())

    plt.figure(figsize=(12, 8))  # Adjusting figure size
    plt.barh(contig_labels, coverage_counts, color='skyblue')
    plt.xlabel('Coverage Count', fontsize=12)
    plt.title('Coverage of Reads Across Merged Contigs', fontsize=14)

    # Add data labels to each bar
    for index, value in enumerate(coverage_counts):
        plt.text(value, index, str(value), va='center', fontsize=10)

    plt.show()


# Check if merged_final_contigs_list is defined before proceeding
if 'merged_final_contigs_list' in locals():
    # Calculate coverage of reads across the merged final contigs
    coverage_dict = calculate_coverage(reads_list, merged_final_contigs_list)

    # Visualize the coverage
    plot_coverage(coverage_dict)
else:
    print("merged_final_contigs_list is not defined, skipping coverage calculation and plotting.")
