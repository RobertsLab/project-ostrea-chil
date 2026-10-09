# Shared settings for the 13-wgbs workflow. Sourced by 13.1–13.3; not run directly.
# All paths are relative to the repository root, like the 11-* scripts.
# Submit jobs from the repo root (or use code/13-wgbs-submit.sh, which does).

# ---- inputs -----------------------------------------------------------------
RAW_DIR="/mmfs1/gscratch/scrubbed/sr320/chil-wgbs"   # 15 x paired 2x151 FASTQs
SAMPLES="data/13-wgbs-samples.tsv"                  # explicit sample -> population table

# Same reference + annotation as the current RNA-seq pipeline (09), so DMRs and
# DE genes share coordinates. Sequence names (Chromosome_1A ...) match the GFF.
GANNET="https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-ostrea-chil/data"
GENOME_DIR="data/13-wgbs-genome"                    # Bismark wants a folder, not a file
GENOME_FA="${GENOME_DIR}/merged_out.fasta"
GFF="${GENOME_DIR}/GN.gene.gff3"

# ---- outputs ----------------------------------------------------------------
OUT="output/13-wgbs"
LOGS="${OUT}/logs"

# ---- parameters -------------------------------------------------------------
# Library is directional (R1 ~1% C, R2 ~6% G), so Bismark runs in default mode.
# 5' trim: 10 bp is a safe default for random-primed/adaptase kits. After the
# first run, check the M-bias plots in ${OUT}/methylation and lower to 0–5 if flat.
TRIM_FRONT1=10
TRIM_FRONT2=10
SCORE_MIN="L,0,-0.6"      # more permissive than Bismark's default L,0,-0.2

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
