#!/bin/bash
# Combine per-chromosome normalize outputs into a single loopQTL object.
# Run after the per-chrom array tasks finish (use --dependency=afterok).
#
#SBATCH --job-name=loopQTL_combine
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/normalize_combine_%j.log

set -euo pipefail

module purge
module load gcc/11.3.0 r/4.4.0

PKG=/project2/rhie_130/suhn/zexunwu/loopQTL
RUN_DIR=/project2/rhie_130/suhn/zexunwu/run_loopQTL    # change me
OUT=${RUN_DIR}/loopQTL_02_normalize_combined.rda

cd $RUN_DIR
Rscript ${PKG}/inst/scripts/combine_normalized.R \
  --pkg $PKG \
  --out $OUT \
  ${RUN_DIR}/chrom_out/normalize_chrom_*.rda
