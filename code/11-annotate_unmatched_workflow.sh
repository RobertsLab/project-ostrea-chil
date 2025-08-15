#!/bin/bash

echo "=== Comprehensive Annotation Workflow for Unmatched Sequences ==="

# Set paths
UNMATCHED_FASTA="output/10-gtf-extraction-blastx/unmatched_sequences/unmatched_transcripts.fasta"
OUTPUT_DIR="output/10-gtf-extraction-blastx/unmatched_sequences"
BLASTDB_DIR="blastdb"
NR_DB="/home/shared/16TB_HDD_01/sam/databases/blastdbs/ncbi-nr-20250429.dmnd"

# Create annotation subdirectories
mkdir -p "${OUTPUT_DIR}/annotation_results"
mkdir -p "${OUTPUT_DIR}/interpro_results"
mkdir -p "${OUTPUT_DIR}/hmmer_results"
mkdir -p "${OUTPUT_DIR}/kegg_results"

echo "Input: ${UNMATCHED_FASTA}"
echo "Output directory: ${OUTPUT_DIR}"

# Check if seqtk is available
if ! command -v seqtk &> /dev/null; then
    echo "Installing seqtk..."
    conda install -c bioconda seqtk -y
fi

echo ""
echo "=== Step 1: Additional BLAST searches ==="

# Check if NR database exists and run BLASTx
if [ -f "${NR_DB}" ]; then
    echo "Running Diamond BLASTx against NR database (${NR_DB})..."
    diamond blastx \
        --query "${UNMATCHED_FASTA}" \
        --db "${NR_DB}" \
        --out "${OUTPUT_DIR}/annotation_results/blastx_nr_results.tsv" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --max-target-seqs 5 \
        --evalue 1e-5 \
        --threads 8 \
        --sensitive
    
    echo "✓ NR BLASTx completed"
else
    echo "NR database not found at ${NR_DB}. Skipping NR search."
fi

# Check if RefSeq database exists
if [ -f "${BLASTDB_DIR}/refseq_protein.dmnd" ]; then
    echo "Running Diamond BLASTx against RefSeq database..."
    diamond blastx \
        --query "${UNMATCHED_FASTA}" \
        --db "${BLASTDB_DIR}/refseq_protein.dmnd" \
        --out "${OUTPUT_DIR}/annotation_results/blastx_refseq_results.tsv" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --max-target-seqs 5 \
        --evalue 1e-5 \
        --threads 8 \
        --sensitive
    
    echo "✓ RefSeq BLASTx completed"
else
    echo "RefSeq database not found. Skipping RefSeq search."
fi

echo ""
echo "=== Step 2: InterProScan Functional Annotation ==="

# Check if InterProScan is available
if command -v interproscan.sh &> /dev/null; then
    echo "Running InterProScan for functional annotation..."
    interproscan.sh \
        -i "${UNMATCHED_FASTA}" \
        -o "${OUTPUT_DIR}/interpro_results" \
        -f tsv \
        -dp \
        --goterms \
        --pathways \
        --cpu 8
    
    echo "✓ InterProScan completed"
else
    echo "InterProScan not found. Install with: conda install -c bioconda interproscan"
fi

echo ""
echo "=== Step 3: HMMER Domain Prediction ==="

# Check if HMMER is available
if command -v hmmsearch &> /dev/null; then
    echo "Running HMMER against Pfam database..."
    
    # Check if Pfam database exists
    if [ -f "${BLASTDB_DIR}/Pfam-A.hmm" ]; then
        hmmsearch \
            --cpu 8 \
            --tblout "${OUTPUT_DIR}/hmmer_results/pfam_results.txt" \
            "${BLASTDB_DIR}/Pfam-A.hmm" \
            "${UNMATCHED_FASTA}"
        
        echo "✓ HMMER Pfam search completed"
    else
        echo "Pfam database not found. Download from: https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
    fi
else
    echo "HMMER not found. Install with: conda install -c bioconda hmmer"
fi

echo ""
echo "=== Step 4: KEGG Pathway Analysis ==="

# Check if KEGG database exists
if [ -f "${BLASTDB_DIR}/kegg_proteins.dmnd" ]; then
    echo "Running Diamond BLASTx against KEGG database..."
    diamond blastx \
        --query "${UNMATCHED_FASTA}" \
        --db "${BLASTDB_DIR}/kegg_proteins.dmnd" \
        --out "${OUTPUT_DIR}/kegg_results/blastx_kegg_results.tsv" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --max-target-seqs 5 \
        --evalue 1e-5 \
        --threads 8 \
        --sensitive
    
    echo "✓ KEGG BLASTx completed"
else
    echo "KEGG database not found. Skipping KEGG search."
fi

echo ""
echo "=== Step 5: Summary and Statistics ==="

# Count sequences in each result file
echo "Annotation Summary:"
echo "=================="

if [ -f "${OUTPUT_DIR}/annotation_results/blastx_nr_results.tsv" ]; then
    NR_HITS=$(cut -f1 "${OUTPUT_DIR}/annotation_results/blastx_nr_results.tsv" | sort -u | wc -l)
    echo "Sequences with NR hits: ${NR_HITS}"
fi

if [ -f "${OUTPUT_DIR}/annotation_results/blastx_refseq_results.tsv" ]; then
    REFSEQ_HITS=$(cut -f1 "${OUTPUT_DIR}/annotation_results/blastx_refseq_results.tsv" | sort -u | wc -l)
    echo "Sequences with RefSeq hits: ${REFSEQ_HITS}"
fi

if [ -f "${OUTPUT_DIR}/interpro_results/*.tsv" ]; then
    INTERPRO_HITS=$(find "${OUTPUT_DIR}/interpro_results" -name "*.tsv" -exec wc -l {} + | tail -1 | awk '{print $1}')
    echo "InterProScan results available"
fi

if [ -f "${OUTPUT_DIR}/hmmer_results/pfam_results.txt" ]; then
    HMMER_HITS=$(grep -v "^#" "${OUTPUT_DIR}/hmmer_results/pfam_results.txt" | wc -l)
    echo "HMMER Pfam hits: ${HMMER_HITS}"
fi

echo ""
echo "=== Next Steps ==="
echo "1. Review annotation results in: ${OUTPUT_DIR}"
echo "2. Combine results from multiple databases"
echo "3. Perform functional enrichment analysis"
echo "4. Create annotation summary report"
echo "5. Visualize results using R/Python"
