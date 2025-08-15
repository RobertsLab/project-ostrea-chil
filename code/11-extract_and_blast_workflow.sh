#!/bin/bash

# Workflow to extract sequences from merged GTF and run diamond blastx
# Author: Generated workflow
# Date: $(date)

set -e  # Exit on any error

echo "=== Starting GTF to FASTA extraction and Diamond BLASTx workflow ==="

# Set paths
GTF_FILE="output/09-align-v1-clean/stringtie_out/merged.gtf"
GENOME_FILE="data/merged_out.fasta"  # Using merged genome assembly
OUTPUT_DIR="output/10-gtf-extraction-blastx"
FASTA_OUTPUT="${OUTPUT_DIR}/merged_transcripts.fasta"
DIAMOND_DB_DIR="${OUTPUT_DIR}/diamond_db"
DIAMOND_OUTPUT="${OUTPUT_DIR}/blastx_results.tsv"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DIAMOND_DB_DIR}"

echo "=== Step 1: Extracting sequences from merged GTF ==="
echo "Input GTF: ${GTF_FILE}"
echo "Reference genome: ${GENOME_FILE}"
echo "Output FASTA: ${FASTA_OUTPUT}"

# Extract sequences using gffread
gffread "${GTF_FILE}" \
    -g "${GENOME_FILE}" \
    -w "${FASTA_OUTPUT}" \
    --force-exons

echo "✓ Sequences extracted successfully"
echo "FASTA file size: $(du -h "${FASTA_OUTPUT}" | cut -f1)"

# Count sequences
SEQUENCE_COUNT=$(grep -c "^>" "${FASTA_OUTPUT}")
echo "Number of sequences extracted: ${SEQUENCE_COUNT}"

echo ""
echo "=== Step 2: Preparing Diamond BLASTx ==="

# Check if Diamond is available
if ! command -v diamond &> /dev/null; then
    echo "❌ Diamond not found. Please install it first:"
    echo "   conda install -c bioconda diamond"
    exit 1
fi

echo "✓ Diamond found: $(diamond version)"

echo ""
echo "=== Step 3: Running Diamond BLASTx against NR database ==="
echo "Note: This step requires the NR database to be downloaded separately."
echo "To download NR database, run:"
echo "   diamond makedb --in nr.gz --db ${DIAMOND_DB_DIR}/nr"
echo ""

# Check if NR database exists
NR_DB="${DIAMOND_DB_DIR}/nr.dmnd"
if [ -f "${NR_DB}" ]; then
    echo "✓ NR database found: ${NR_DB}"
    
    echo "Running Diamond BLASTx..."
    diamond blastx \
        --query "${FASTA_OUTPUT}" \
        --db "${NR_DB}" \
        --out "${DIAMOND_OUTPUT}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle \
        --max-target-seqs 5 \
        --evalue 1e-5 \
        --threads 8 \
        --sensitive
    
    echo "✓ Diamond BLASTx completed successfully"
    echo "Results saved to: ${DIAMOND_OUTPUT}"
    
    # Show summary statistics
    echo ""
    echo "=== BLASTx Results Summary ==="
    echo "Total queries: ${SEQUENCE_COUNT}"
    echo "Total hits: $(wc -l < "${DIAMOND_OUTPUT}")"
    echo "Unique queries with hits: $(cut -f1 "${DIAMOND_OUTPUT}" | sort -u | wc -l)"
    
else
    echo "❌ NR database not found at: ${NR_DB}"
    echo ""
    echo "To download and prepare the NR database:"
    echo "1. Download NR database (this may take a while):"
    echo "   wget -O ${DIAMOND_DB_DIR}/nr.gz https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
    echo ""
    echo "2. Create Diamond database:"
    echo "   diamond makedb --in ${DIAMOND_DB_DIR}/nr.gz --db ${NR_DB}"
    echo ""
    echo "3. Then re-run this script to perform BLASTx"
fi

echo ""
echo "=== Workflow completed ==="
echo "Output files:"
echo "  - Extracted sequences: ${FASTA_OUTPUT}"
echo "  - Diamond database directory: ${DIAMOND_DB_DIR}"
if [ -f "${DIAMOND_OUTPUT}" ]; then
    echo "  - BLASTx results: ${DIAMOND_OUTPUT}"
fi
