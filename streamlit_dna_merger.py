import streamlit as st
from Bio import Seq
from Bio import SeqIO
from difflib import ndiff
import io

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    return str(Seq.Seq(seq).reverse_complement())

def find_overlap(seq1, seq2, min_overlap=10):
    """Finds the maximum overlap between two sequences and returns the merged sequence."""
    max_overlap = 0
    merged_sequence = ""
    best_match = ""
    
    for i in range(min_overlap, len(seq1)):
        if seq2.startswith(seq1[i:]):
            max_overlap = len(seq1) - i
            merged_sequence = seq1 + seq2[max_overlap:]
            best_match = seq1[i:]
            break
    
    return merged_sequence, max_overlap, best_match

def identify_differences(seq1, seq2):
    """Identifies differences in the overlapping region."""
    differences = []
    for i, s in enumerate(ndiff(seq1, seq2)):
        if s[0] != ' ':  # Indicates a difference
            differences.append(f"Position {i}: {s}")
    return differences

def merge_sequences(sequences):
    """Merges a list of DNA sequences, considering sense and antisense directions."""
    merged_seq = sequences[0]
    results = []
    for i in range(1, len(sequences)):
        sense_seq = sequences[i]
        antisense_seq = reverse_complement(sequences[i])
        
        merged_sense, overlap_sense, match_sense = find_overlap(merged_seq, sense_seq)
        merged_antisense, overlap_antisense, match_antisense = find_overlap(merged_seq, antisense_seq)
        
        if overlap_sense >= overlap_antisense:
            merged_seq = merged_sense
            best_match = match_sense
            overlap = overlap_sense
            orientation = "Sense"
        else:
            merged_seq = merged_antisense
            best_match = match_antisense
            overlap = overlap_antisense
            orientation = "Antisense"
        
        percent_identity = (overlap / len(best_match)) * 100 if best_match else 0
        differences = identify_differences(best_match, sequences[i][:overlap])
        
        results.append((i, overlap, percent_identity, orientation, differences))
    
    return merged_seq, results

def parse_fasta(fasta_text):
    """Parses FASTA format sequences from user input."""
    sequences = []
    fasta_io = io.StringIO(fasta_text)  # Convert string to file-like object
    for record in SeqIO.parse(fasta_io, "fasta"):
        sequences.append(str(record.seq).upper())
    return sequences

st.title("DNA Sequence Merger")

fasta_input = st.text_area("Paste FASTA sequences here:")

if st.button("Merge Sequences"):
    if fasta_input.strip():
        sequences = parse_fasta(fasta_input)
        if len(sequences) < 2:
            st.error("Please provide at least two sequences in FASTA format.")
        else:
            merged_sequence, merge_results = merge_sequences(sequences)
            st.subheader("Merged Sequence:")
            st.text_area("", merged_sequence, height=200)
            
            st.subheader("Merge Details:")
            for res in merge_results:
                i, overlap, percent_identity, orientation, differences = res
                st.write(f"**Merged with sequence {i}:**")
                st.write(f"Overlap: {overlap} bases, Percent Identity: {percent_identity:.2f}%, Orientation: {orientation}")
                if differences:
                    st.write("Differences in overlap:")
                    st.write("\n".join(differences))
    else:
        st.error("Please paste sequences in FASTA format.")
