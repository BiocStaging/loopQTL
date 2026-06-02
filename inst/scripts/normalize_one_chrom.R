#!/usr/bin/env Rscript
#
# Per-chromosome loopQTL_normalize() driver, designed for a SLURM job array.
#
# Each array task runs normalize for exactly ONE chromosome and saves a
# per-chrom .rda. A final combine job loads them all and stitches them.
#
# NOTE: As of the RINT-moves-to-run_QTL refactor, loopQTL_normalize() no
# longer applies RINT. The phenotype slot stores pre-RINT BNBC/fastlo
# output. Apply RINT (or not) via `normal = TRUE/FALSE` in run_QTL().
#
# Usage (inside a SLURM array task):
#   Rscript normalize_one_chrom.R \
#     --setup    loopQTL_01_setup.rda \
#     --chrom    "$CHROM" \
#     --method   multiHiCcompare \
#     --span     0.7 \
#     --pkg      /project2/rhie_130/suhn/zexunwu/loopQTL \
#     --out      normalize_chrom_${SLURM_ARRAY_TASK_ID}.rda
#
# CHROM can come from either:
#   - $SLURM_ARRAY_TASK_ID indexing into a chrom list you build below, OR
#   - a CHROMS file you sed by $SLURM_ARRAY_TASK_ID.

suppressPackageStartupMessages({
  library(BiocParallel)
})

# ---- minimal arg parser (no extra deps) ---------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) return(default)
  args[i + 1L]
}

setup_path <- get_arg("setup")
chrom      <- get_arg("chrom")
method     <- get_arg("method",  "multiHiCcompare")
span       <- as.numeric(get_arg("span", "0.7"))
pkg_path   <- get_arg("pkg")
out_path   <- get_arg("out")

stopifnot(
  !is.null(setup_path), !is.null(chrom),
  !is.null(pkg_path),   !is.null(out_path),
  file.exists(setup_path),
  dir.exists(pkg_path)
)

devtools::load_all(pkg_path)

message("[", Sys.time(), "] loading setup obj from ", setup_path)
load(setup_path)   # restores `obj`
stopifnot(exists("obj"))

message("[", Sys.time(), "] normalizing chrom ", chrom,
        " (method=", method, ", span=", span, ")")
obj <- loopQTL_normalize(
  obj,
  method  = method,
  span    = span,
  chroms  = chrom,
  BPPARAM = BiocParallel::SerialParam()   # one chrom, one process
)

n_filled <- sum(rowSums(!is.na(obj@phenotype)) > 0L)
message("[", Sys.time(), "] DONE; ", n_filled, " loop rows filled. Saving to ", out_path)
save(obj, file = out_path)
