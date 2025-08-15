#!/bin/bash

echo "=== Extracting sequences from merged GTF ==="

# Set paths
GTF_FILE="output/09-align-v1-clean/stringtie_out/merged.gtf"
GENOME_FILE="data/merged_out.fasta"
OUTPUT_DIR="output/10-gtf-extraction-blastx"
FASTA_OUTPUT="${OUTPUT_DIR}/merged_transcripts.fasta"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Input GTF: ${GTF_FILE}"
echo "Reference genome: ${GENOME_FILE}"
echo "Output FASTA: ${FASTA_OUTPUT}"

# Extract sequences using gffread
echo "Running gffread..."
gffread "${GTF_FILE}" \
    -g "${GENOME_FILE}" \
    -w "${FASTA_OUTPUT}" \
    --force-exons

echo "✓ Sequences extracted successfully!"
echo "FASTA file size: $(du -h "${FASTA_OUTPUT}" | cut -f1)"

# Count sequences
SEQUENCE_COUNT=$(grep -c "^>" "${FASTA_OUTPUT}")
echo "Number of sequences extracted: ${SEQUENCE_COUNT}"

# Show first few sequences
echo ""
echo "First few sequences:"
head -20 "${FASTA_OUTPUT}"

echo ""
echo "=== Next steps ==="
echo "Swiss-Prot Diamond database found in blastdb/ directory"
echo "Running Diamond BLASTx against Swiss-Prot database..."
echo ""

# Run Diamond BLASTx against Swiss-Prot database
diamond blastx \
    --query "${FASTA_OUTPUT}" \
    --db "blastdb/uniprot_sprot_r2024_03_v0.9.19.dmnd" \
    --out "${OUTPUT_DIR}/blastx_swissprot_results.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --max-target-seqs 5 \
    --evalue 1e-5 \
    --threads 8 \
    --sensitive

echo ""
echo "✓ Diamond BLASTx completed successfully!"
echo "Results saved to: ${OUTPUT_DIR}/blastx_swissprot_results.tsv"

# Show summary statistics
echo ""
echo "=== BLASTx Results Summary ==="
echo "Total queries: ${SEQUENCE_COUNT}"
echo "Total hits: $(wc -l < "${OUTPUT_DIR}/blastx_swissprot_results.tsv")"
echo "Unique queries with hits: $(cut -f1 "${OUTPUT_DIR}/blastx_swissprot_results.tsv" | sort -u | wc -l)"
