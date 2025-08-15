#!/bin/bash

echo "=== Simple Annotation Workflow for Unmatched Sequences ==="

# Set paths
UNMATCHED_FASTA="output/10-gtf-extraction-blastx/unmatched_sequences/unmatched_transcripts.fasta"
OUTPUT_DIR="output/10-gtf-extraction-blastx/unmatched_sequences"
SWISSPROT_DB="blastdb/uniprot_sprot_r2024_03_v0.9.19.dmnd"

# Create output directories
mkdir -p "${OUTPUT_DIR}/annotation_results"
mkdir -p "${OUTPUT_DIR}/sequence_analysis"

echo "Input: ${UNMATCHED_FASTA}"
echo "Output directory: ${OUTPUT_DIR}"

# Check sequence statistics
echo ""
echo "=== Sequence Analysis ==="
TOTAL_SEQUENCES=$(grep -c "^>" "${UNMATCHED_FASTA}")
echo "Total unmatched sequences: ${TOTAL_SEQUENCES}"

# Analyze sequence lengths
echo "Analyzing sequence lengths..."
awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq = seq $0} END {if (seq) print length(seq)}' "${UNMATCHED_FASTA}" > "${OUTPUT_DIR}/sequence_analysis/sequence_lengths.txt"

# Calculate length statistics
echo "Sequence length statistics:"
echo "Min length: $(sort -n "${OUTPUT_DIR}/sequence_analysis/sequence_lengths.txt" | head -1)"
echo "Max length: $(sort -n "${OUTPUT_DIR}/sequence_analysis/sequence_lengths.txt" | tail -1)"
echo "Average length: $(awk '{sum+=$1} END {print sum/NR}' "${OUTPUT_DIR}/sequence_analysis/sequence_lengths.txt" | awk '{printf "%.1f", $1}')"

# Check for potential ORFs
echo ""
echo "=== ORF Analysis ==="
echo "Checking for potential open reading frames..."

# Simple ORF detection (basic approach)
awk '/^>/ {if (seq) {
    # Count ATG, TAA, TAG, TGA
    atg_count = gsub(/ATG/, "ATG", seq)
    taa_count = gsub(/TAA/, "TAA", seq)
    tag_count = gsub(/TAG/, "TAG", seq)
    tga_count = gsub(/TGA/, "TGA", seq)
    total_stops = taa_count + tag_count + tga_count
    
    print $0 "\t" atg_count "\t" total_stops "\t" length(seq)
}
seq=""; next} {seq = seq $0} END {
    if (seq) {
        atg_count = gsub(/ATG/, "ATG", seq)
        taa_count = gsub(/TAA/, "TAA", seq)
        tag_count = gsub(/TAG/, "TAG", seq)
        tga_count = gsub(/TGA/, "TGA", seq)
        total_stops = taa_count + tag_count + tga_count
        print $0 "\t" atg_count "\t" total_stops "\t" length(seq)
    }
}' "${UNMATCHED_FASTA}" > "${OUTPUT_DIR}/sequence_analysis/orf_analysis.txt"

echo "ORF analysis saved to: ${OUTPUT_DIR}/sequence_analysis/orf_analysis.txt"

# Try to run Diamond with relaxed parameters against Swiss-Prot
echo ""
echo "=== Re-running Swiss-Prot Search with Relaxed Parameters ==="
echo "Running Diamond BLASTx with more relaxed E-value threshold..."

diamond blastx \
    --query "${UNMATCHED_FASTA}" \
    --db "${SWISSPROT_DB}" \
    --out "${OUTPUT_DIR}/annotation_results/blastx_swissprot_relaxed.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --max-target-seqs 10 \
    --evalue 1e-3 \
    --threads 8 \
    --sensitive

echo "✓ Relaxed Swiss-Prot search completed"

# Count additional hits
RELAXED_HITS=$(cut -f1 "${OUTPUT_DIR}/annotation_results/blastx_swissprot_relaxed.tsv" | sort -u | wc -l)
echo "Sequences with relaxed Swiss-Prot hits: ${RELAXED_HITS}"

# Create summary report
echo ""
echo "=== Annotation Summary Report ==="
echo "================================="
echo "Total sequences analyzed: ${TOTAL_SEQUENCES}"
echo "Original Swiss-Prot hits: 17,026"
echo "Relaxed Swiss-Prot hits: ${RELAXED_HITS}"
echo "Still unmatched: $((TOTAL_SEQUENCES - RELAXED_HITS))"

# Save summary to file
cat > "${OUTPUT_DIR}/annotation_summary.txt" << EOF
Annotation Summary Report
=========================
Date: $(date)
Total sequences analyzed: ${TOTAL_SEQUENCES}
Original Swiss-Prot hits: 17,026
Relaxed Swiss-Prot hits: ${RELAXED_HITS}
Still unmatched: $((TOTAL_SEQUENCES - RELAXED_HITS))

Recommendations for remaining unmatched sequences:
1. Run against NCBI NR database (requires compatible Diamond version)
2. Use InterProScan for domain prediction
3. Perform de novo functional prediction
4. Analyze expression patterns and co-expression networks
5. Compare with related species genomes
EOF

echo ""
echo "=== Summary saved to: ${OUTPUT_DIR}/annotation_summary.txt ==="
echo ""
echo "=== Next Steps ==="
echo "1. Review the annotation summary"
echo "2. Consider installing InterProScan for domain prediction"
echo "3. Analyze expression patterns of unmatched sequences"
echo "4. Compare with related species for evolutionary insights"
echo "5. Use machine learning approaches for functional prediction"
