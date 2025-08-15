#!/bin/bash

echo "=== Identifying sequences without Swiss-Prot matches ==="

# Set paths
FASTA_FILE="output/10-gtf-extraction-blastx/merged_transcripts.fasta"
BLASTX_RESULTS="output/10-gtf-extraction-blastx/blastx_swissprot_results.tsv"
OUTPUT_DIR="output/10-gtf-extraction-blastx"

# Create output directory for unmatched sequences
mkdir -p "${OUTPUT_DIR}/unmatched_sequences"

echo "Input FASTA: ${FASTA_FILE}"
echo "BLASTx results: ${BLASTX_RESULTS}"

# Get list of sequences that had hits
echo "Extracting sequences with Swiss-Prot hits..."
cut -f1 "${BLASTX_RESULTS}" | sort -u > "${OUTPUT_DIR}/matched_sequences.txt"

# Get total sequence count
TOTAL_SEQUENCES=$(grep -c "^>" "${FASTA_FILE}")
MATCHED_SEQUENCES=$(wc -l < "${OUTPUT_DIR}/matched_sequences.txt")
UNMATCHED_COUNT=$((TOTAL_SEQUENCES - MATCHED_SEQUENCES))

echo "Total sequences: ${TOTAL_SEQUENCES}"
echo "Sequences with Swiss-Prot hits: ${MATCHED_SEQUENCES}"
echo "Sequences without Swiss-Prot hits: ${UNMATCHED_COUNT}"

# Extract unmatched sequences
echo "Extracting unmatched sequences..."
grep "^>" "${FASTA_FILE}" | sed 's/^>//' > "${OUTPUT_DIR}/all_sequence_ids.txt"

# Find unmatched sequences
comm -23 <(sort "${OUTPUT_DIR}/all_sequence_ids.txt") <(sort "${OUTPUT_DIR}/matched_sequences.txt") > "${OUTPUT_DIR}/unmatched_sequence_ids.txt"

# Extract unmatched sequences to FASTA
echo "Creating FASTA file of unmatched sequences..."
seqtk subseq "${FASTA_FILE}" "${OUTPUT_DIR}/unmatched_sequence_ids.txt" > "${OUTPUT_DIR}/unmatched_sequences/unmatched_transcripts.fasta"

echo "✓ Unmatched sequences extracted successfully!"
echo "Unmatched sequences saved to: ${OUTPUT_DIR}/unmatched_sequences/unmatched_transcripts.fasta"
echo ""

echo "=== Next annotation steps for unmatched sequences ==="
echo "1. Run against NCBI NR database (if available)"
echo "2. Run against RefSeq database"
echo "3. Run against custom protein databases"
echo "4. Functional annotation using InterProScan"
echo "5. Domain prediction using HMMER"
echo "6. Gene ontology prediction"
echo "7. Pathway analysis using KEGG"
echo "8. De novo functional prediction using machine learning"
