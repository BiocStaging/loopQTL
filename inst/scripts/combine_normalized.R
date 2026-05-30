#!/usr/bin/env Rscript
#
# Combine per-chromosome normalize outputs into a single loopQTL object.
#
# Usage:
#   Rscript combine_normalized.R \
#     --pkg /project2/rhie_130/suhn/zexunwu/loopQTL \
#     --out loopQTL_02_normalize_combined.rda \
#     normalize_chrom_*.rda

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) return(default)
  args[i + 1L]
}

pkg_path <- get_arg("pkg")
out_path <- get_arg("out")

# Anything that's not a --flag and not the value after a --flag is treated as a path
flag_positions <- grep("^--", args)
value_positions <- flag_positions + 1L
non_path_idx <- sort(unique(c(flag_positions, value_positions)))
paths <- args[setdiff(seq_along(args), non_path_idx)]

stopifnot(
  !is.null(pkg_path), !is.null(out_path), length(paths) > 0L,
  dir.exists(pkg_path)
)

devtools::load_all(pkg_path)

message("[", Sys.time(), "] combining ", length(paths), " per-chrom objects")
obj <- loopQTL_combine_normalized(paths = paths)
message("[", Sys.time(), "] chroms_normalized: ",
        paste(obj@meta$chroms_normalized, collapse = ", "))
message("[", Sys.time(), "] saving to ", out_path)
save(obj, file = out_path)
