#!/bin/bash
# Per-chromosome loopQTL_normalize as a SLURM array.
#
# Each array task processes ONE chromosome.
# Each task gets the FULL per-node memory budget for just that chrom.
#
# Submit with:
#   sbatch --array=1-23 normalize_array.slurm
# (1..22 for autosomes; add chrX as task 23.)
#
# After the array finishes, kick off the combine job:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> normalize_combine.slurm
#
#SBATCH --job-name=loopQTL_norm
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=logs/norm_chrom_%A_%a.log

set -euo pipefail

module purge
module load gcc/11.3.0 r/4.4.0

PKG=/project2/rhie_130/suhn/zexunwu/loopQTL
RUN_DIR=/project2/rhie_130/suhn/zexunwu/run_loopQTL    # change me
SETUP=${RUN_DIR}/loopQTL_01_setup.rda

# Map SLURM_ARRAY_TASK_ID -> chrom name. Adjust naming if your .hic uses "chr1" vs "1".
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHROM=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}

mkdir -p ${RUN_DIR}/logs ${RUN_DIR}/chrom_out
OUT=${RUN_DIR}/chrom_out/normalize_chrom_${CHROM}.rda

cd $RUN_DIR
Rscript ${PKG}/inst/scripts/normalize_one_chrom.R \
  --setup  $SETUP \
  --chrom  $CHROM \
  --method BNBC \
  --pkg    $PKG \
  --out    $OUT
