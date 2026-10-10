# Shared settings for the 13-wgbs workflow. Sourced by 13.1–13.3; not run directly.
# All paths are relative to the repository root, like the 11-* scripts.
# Submit jobs from the repo root (or use code/13-wgbs-submit.sh, which does).

# ---- inputs -----------------------------------------------------------------
RAW_DIR="/mmfs1/gscratch/scrubbed/sr320/chil-wgbs"   # 15 x paired 2x151 FASTQs
SAMPLES="data/13-wgbs-samples.tsv"                  # explicit sample -> population table

# Reference, chosen with WGBS_REF (default "merged"):
#   merged  merged_out.fasta + GN.gene.gff3, the same reference as the RNA-seq
#           (09), so DMRs and DE genes share coordinates. It is HapA's
#           chromosomes 2A/4A/6A/10A plus HapB's 1B/5B/6B/7B/8B/9B.
#   hapa    Och_HapA_assembly.fa alone (no annotation), to test whether results
#           depend on the mixed reference. Separate genome folder and outputs:
#           WGBS_REF=hapa bash code/13-wgbs-submit.sh
# GENOME_FILES lists "url<TAB>filename" pairs that 13.1 downloads.
GANNET="https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-ostrea-chil/data"
WGBS_REF="${WGBS_REF:-merged}"
case "${WGBS_REF}" in
  merged)
    GENOME_DIR="data/13-wgbs-genome"                # Bismark wants a folder, not a file
    GENOME_FA="${GENOME_DIR}/merged_out.fasta"
    GFF="${GENOME_DIR}/GN.gene.gff3"
    GENOME_FILES="${GANNET}/merged_out.fasta	merged_out.fasta
${GANNET}/GN.gene.gff3	GN.gene.gff3"
    OUT="output/13-wgbs"
    ;;
  hapa)
    GENOME_DIR="data/13-wgbs-genome-hapa"
    GENOME_FA="${GENOME_DIR}/Och_HapA_assembly.fa"
    GFF=""                                          # no annotation on HapA coordinates
    GENOME_FILES="https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/OCEAN/docs/jbrowse/data/HapA/Och_HapA_assembly.fa	Och_HapA_assembly.fa"
    OUT="output/13-wgbs-hapa"
    ;;
  *)
    echo "ERROR: WGBS_REF must be merged or hapa, not ${WGBS_REF}" >&2
    exit 1
    ;;
esac

# ---- outputs ----------------------------------------------------------------
LOGS="${OUT}/logs"

# ---- parameters -------------------------------------------------------------
# Library is directional (R1 ~1% C, R2 ~6% G), so Bismark runs in default mode.
# 5' trim: 10 bp is a safe default for random-primed/adaptase kits. After the
# first run, check the M-bias plots in ${OUT}/methylation and lower to 0–5 if flat.
TRIM_FRONT1=10
TRIM_FRONT2=10
SCORE_MIN="L,0,-0.6"      # more permissive than Bismark's default L,0,-0.2
# Pum1 test M-bias: R1 flat after the 10 bp trim, but R2 still elevated at its
# first 5 bases (29% -> 20% CpG meth, flat from base 6). Skip those 5 bases at
# extraction instead of re-trimming, so no re-alignment is needed.
IGNORE_R2=5

# ---- SLURM ----------------------------------------------------------------------
# Applied by 13-wgbs-submit.sh (overrides the #SBATCH lines in each script).
# coenv on cpu-g2 (192-core / 1.5 TB nodes) is a shared 928-core allocation.
# Kept polite: 5 alignment tasks at a time x 32 cores = 160 cores (~17% of coenv).
# Lower if coenv is busy; or switch back to srlab / cpu-g2-mem2x.
SLURM_ACCOUNT="coenv"
SLURM_PARTITION="cpu-g2"
ALIGN_CONCURRENCY=5

# ---- software -----------------------------------------------------------------
# The srlab `myflow` env has fastp 1.0.1, Bismark 0.25.1, bowtie2, samtools 1.22,
# multiqc. If you move clusters, this is the only line to change.
CONDA_ROOT="/mmfs1/gscratch/srlab/sr320/miniforge3"
activate_env() {
  # shellcheck disable=SC1091
  source "${CONDA_ROOT}/etc/profile.d/conda.sh"
  conda activate myflow
  # The env was moved from ~/miniforge3, but its curl still has the old
  # CA-bundle path compiled in; point it at the env's own certificates.
  export CURL_CA_BUNDLE="${CONDA_PREFIX}/ssl/cacert.pem"
}

# Refuse to run from anywhere but the repo root, so relative paths are right.
if [[ ! -f "${SAMPLES}" ]]; then
  echo "ERROR: ${SAMPLES} not found - run from the repository root." >&2
  exit 1
fi
