import numpy as np
from Bio import SeqIO

# Motif matrix for scoring sequences
motif = np.array([
    [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 2, -99, -99, 0.5],  # Scores for A
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],  # Scores for T
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],  # Scores for C
    [0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0.5, -99, 2, 0]  # Scores for G
])

# Assign indices for each base
base_of_index = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
min_ORF_len = 60  # Minimum length of ORF
cut_off = 7.25  # Score cutoff for identifying ORF


# scoring the first 13 bases of a sequence
def scoreMotif(sequence):
    motifscore = 0
    for i, base in enumerate(sequence[:13]):
        if base in base_of_index:  # Check if base is valid
            motifscore += motif[base_of_index[base], i]
    return motifscore


# scaning ORFs in a sequence
def scanSeq(sequence):
    potentialStartPos = []  # Stores ORF start positions
    ORFlengths = []  # Stores ORF lengths
    ORFSeqs = []  # Stores ORF sequences

    location = set()  # To avoid duplicate detection

    # finding start codon
    def orf(id):
        y = -np.inf
        z = None
        for i in range(id, id + 13, 3):
            codon = sequence[i:i + 3]
            if codon in ('ATG', 'GTG'):
                motifscore = scoreMotif(sequence[i:i + 13])
                if motifscore > y:
                    y = motifscore
                    z = i
        return z

    # Scaning the whole sequence
    for i in range(len(sequence) - 12):
        if scoreMotif(sequence[i:i + 13]) > cut_off:
            b = orf(i)
            if b is not None and b not in location:
                location.add(b)
                for end in range(b + 3, len(sequence), 3):
                    if sequence[end:end + 3] in ('TAA', 'TAG', 'TGA'):  # Stop codons
                        if end - b + 3 >= min_ORF_len:
                            potentialStartPos.append(b)
                            ORFlengths.append(end - b + 3)
                            ORFSeqs.append(sequence[b:end + 3])
                        break
    return potentialStartPos, ORFlengths, ORFSeqs


# reading the FASTA file using SeqIO
def reading_fasta(filepath):
    contigs = {}
    for record in SeqIO.parse(filepath, "fasta"):
        contigs[record.id] = str(record.seq)  # Store ID and sequence
    return contigs


# identify ORFs in a set of contigs
def identifyORFs(contigs):
    finalseq = []
    for contig, seq in contigs.items():
        st, lengths, read = scanSeq(seq)
        for t, l, s in zip(st, lengths, read):
            finalseq.append(f"> {contig}|Length {l}|at position {t + 1}\n{s}")
    return finalseq


# Main function to run the script
if __name__ == "__main__":
    contig_data = reading_fasta("spaceSeq.fa")  # Read the FASTA file
    orfResults = identifyORFs(contig_data)  # Identify ORFs

    # Write the results to an output FASTA file
    with open("rima.fa", "w") as output_file:
        for result in orfResults:
            output_file.write(result + "\n")

    print("ORF found! Result file: 'rima.fa'.")
