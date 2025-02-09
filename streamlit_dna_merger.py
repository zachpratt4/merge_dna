import streamlit as st
from Bio import Seq
from Bio import SeqIO
from difflib import ndiff
import io

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    return str(Seq.Seq(seq).reverse_complement())

def find_overlap(seq1, seq2, min_overlap=20):
    """Finds the maximum overlap between two sequences with no mismatches allowed."""
    max_overlap = 0
    merged_sequence = ""
    best_match = ""
    
    for i in range(min_overlap, len(seq1)):
        overlap_seq1 = seq1[i:]
        overlap_seq2 = seq2[:len(overlap_seq1)]
        
        if overlap_seq1 == overlap_seq2:
            max_overlap = len(overlap_seq1)
            merged_sequence = seq1 + seq2[max_overlap:]
            best_match = overlap_seq1
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
    merged_seq = sequences[0]["sequence"]
    results = []
    for i in range(1, len(sequences)):
        sequence_entry = sequences[i]
        sequence = sequence_entry["sequence"]
        direction = sequence_entry["direction"]
        
        if direction == "antisense":
            sequence = reverse_complement(sequence)
        
        merged_sense, overlap_sense, match_sense = find_overlap(merged_seq, sequence)
        
        if overlap_sense < 20:
            results.append((i, "Overlap too short (<20 bases)", 0, "N/A", []))
            continue
        
        merged_seq = merged_sense
        best_match = match_sense
        overlap = overlap_sense
        orientation = "Sense" if direction == "sense" else "Antisense"
        
        percent_identity = (overlap / len(best_match)) * 100 if best_match else 0
        differences = identify_differences(best_match, sequence[:overlap])
        
        results.append((i, overlap, percent_identity, orientation, differences))
    
    return merged_seq, results

def parse_fasta(fasta_text):
    """Parses FASTA format sequences from user input with direction information."""
    sequences = []
    fasta_io = io.StringIO(fasta_text)  # Convert string to file-like object
    for record in SeqIO.parse(fasta_io, "fasta"):
        header = record.description.split()
        direction = "sense"
        if len(header) > 1 and header[1].lower() == "antisense":
            direction = "antisense"
        sequences.append({"sequence": str(record.seq).upper(), "direction": direction})
    return sequences

st.title("DNA Sequence Merger")

st.write("""
Paste FASTA sequences below. If a sequence is in the antisense direction, specify it in the header by adding 'antisense' after the sequence name.
Example:
```
>sequence1
ATGCGTACGTTAGC
>sequence2 antisense
GCTAACGTACGCAT
```
The merged DNA sequence will be output in the 5' to 3' orientation, along with its length.
""")

fasta_input = st.text_area("Paste FASTA sequences here:")

if st.button("Merge Sequences"):
    if fasta_input.strip():
        sequences = parse_fasta(fasta_input)
        if len(sequences) < 2:
            st.error("Please provide at least two sequences in FASTA format.")
        else:
            merged_sequence, merge_results = merge_sequences(sequences)
            if merged_sequence is None:
                st.error("Some sequences had an overlap too short (<20 bases). Unable to merge all sequences.")
                for res in merge_results:
                    i, overlap, _, _, _ = res
                    if isinstance(overlap, str):
                        st.write(f"**Sequence {i}: {overlap}**")
            else:
                st.subheader("Merged Sequence:")
                st.text_area("", merged_sequence, height=200)
                st.write(f"**Length: {len(merged_sequence)} bases**")
                
                st.subheader("Merge Details:")
                for res in merge_results:
                    i, overlap, percent_identity, orientation, differences = res
                    if isinstance(overlap, str):
                        st.write(f"**Sequence {i}: {overlap}**")
                    else:
                        st.write(f"**Merged with sequence {i}:**")
                        st.write(f"Overlap: {overlap} bases, Percent Identity: {percent_identity:.2f}%, Orientation: {orientation}")
                        if differences:
                            st.write("Differences in overlap:")
                            st.write("\n".join(differences))
    else:
        st.error("Please paste sequences in FASTA format.")
