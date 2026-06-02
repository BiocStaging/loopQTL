#' Normalize loopQTL contacts at whole-genome scale
#'
#' For each chromosome, pulls whole-genome `(x, y, counts)` from every .hic via
#' `strawr::straw("NONE", ...)`, then normalizes across samples:
#'
#'  * `method = "multiHiCcompare"`: builds a per-chrom `hicexp` (groups all = 1)
#'    and runs `multiHiCcompare::fastlo(...)`. `fastlo` is multiHiCcompare's
#'    fast loess normalization, distinct from `cyclic_loess`. We do not use
#'    the differential side of multiHiCcompare.
#'  * `method = "BNBC"`: per-chrom log-transformed contact triplets,
#'    band-wise quantile-normalized across samples via
#'    `preprocessCore::normalize.quantiles()` (Bolstad et al. 2003), as in BNBC
#'    (Fletez-Brant et al. 2021). ComBat batch correction is intentionally
#'    skipped. The sparse triplet data structure is our own implementation;
#'    the `bnbc` Bioconductor package is not used (its dense internals don't
#'    scale to 50+ samples at 10 kb).
#'
#' Normalized values are indexed at the consensus loop bin-pairs to build the
#' (loops x samples) phenotype matrix. **No RINT is applied here** — the
#' phenotype slot stores the BNBC/fastlo output directly. RINT moved to
#' `run_QTL()` (`normal` argument there) so the un-RINT'd matrix remains
#' inspectable for QC.
#'
#' Per-chromosome work is dispatched via `BiocParallel`.
#'
#' @param obj A [loopQTL-class] from [loopQTL_setup()].
#' @param method Normalization method. Default `"multiHiCcompare"`.
#' @param depth_correct For `method = "BNBC"`: apply logCPM depth normalization
#'   before band-wise QN, following Fletez-Brant et al. 2021. Default `TRUE`.
#'   Requires `obj@cis` to be populated (either via the samplesheet at setup
#'   or by setup's auto-compute). Ignored for `method = "multiHiCcompare"`
#'   (fastlo handles depth implicitly).
#' @param span Loess smoothing span passed to `multiHiCcompare::fastlo()`.
#'   Default 0.7 (matches `fastlo`'s default).
#' @param chroms Optional character vector of chromosome names to normalize.
#'   Default `NULL` = all chromosomes present in `obj@loops`. When a subset is
#'   passed, only those chroms' loop rows are filled in `obj@phenotype`; rows
#'   for other chroms stay `NA`. Use this to split a large run across SLURM
#'   array tasks (one chrom per task), then combine the per-chrom objects with
#'   [loopQTL_combine_normalized()].
#' @param BPPARAM A [BiocParallel::BiocParallelParam] for parallelization over
#'   chromosomes. Default [BiocParallel::SerialParam()].
#'
#' @return The input object with `phenotype` filled (BNBC or fastlo output,
#'   pre-RINT) and `meta$normalize_method` recorded.
#'
#' @export
loopQTL_normalize <- function(obj,
                              method        = c("multiHiCcompare", "BNBC"),
                              depth_correct = TRUE,
                              span          = 0.7,
                              chroms        = NULL,
                              BPPARAM       = BiocParallel::SerialParam()) {
  stopifnot(methods::is(obj, "loopQTL"))
  method <- match.arg(method)

  # ---- Depth-correction requires cis when using BNBC ----
  if (method == "BNBC" && isTRUE(depth_correct)) {
    if (length(obj@cis) == 0L) {
      stop("BNBC with depth_correct = TRUE requires per-sample cis depth in ",
           "obj@cis, but it is empty. Options:\n",
           "  1. Include a `cis` column in the samplesheet at loopQTL_setup(), or\n",
           "  2. Re-run loopQTL_setup() with `compute_cis = TRUE` (default), or\n",
           "  3. Manually attach: obj@cis <- setNames(<numeric>, obj@samples), or\n",
           "  4. Set depth_correct = FALSE (not recommended — known to leave\n",
           "     residual depth bias in 50+ sample cohorts).")
    }
    if (length(obj@cis) != length(obj@samples) ||
        !identical(names(obj@cis), obj@samples)) {
      stop("obj@cis must be a named numeric vector with one entry per sample, ",
           "names matching obj@samples")
    }
  }

  samples_vec <- obj@samples
  resolution  <- obj@resolution
  loops_gi    <- obj@loops

  # Per-sample per-chrom hic paths from meta (we re-derive from setup invariants)
  if (is.null(obj@meta$hic_paths)) {
    stop("obj@meta$hic_paths missing; loopQTL_setup must record the per-sample .hic paths")
  }
  hic_paths <- obj@meta$hic_paths

  # Consensus loops keyed by chrom + canonical (binA <= binB)
  anchorA <- InteractionSet::anchors(loops_gi, type = "first")
  anchorB <- InteractionSet::anchors(loops_gi, type = "second")
  loop_chr  <- as.character(GenomicRanges::seqnames(anchorA))
  loop_binA <- GenomicRanges::start(anchorA) - 1L   # 0-based bin start, matches straw
  loop_binB <- GenomicRanges::start(anchorB) - 1L
  swap <- loop_binA > loop_binB
  if (any(swap)) {
    tmp <- loop_binA[swap]; loop_binA[swap] <- loop_binB[swap]; loop_binB[swap] <- tmp
  }
  loop_key <- paste(loop_chr, loop_binA, loop_binB, sep = ":")

  all_chroms <- unique(loop_chr)
  if (is.null(chroms)) {
    chroms_to_run <- all_chroms
  } else {
    chroms_to_run <- intersect(as.character(chroms), all_chroms)
    missing <- setdiff(as.character(chroms), all_chroms)
    if (length(missing) > 0L)
      warning("These `chroms` are not present in obj@loops and will be ignored: ",
              paste(missing, collapse = ", "))
    if (length(chroms_to_run) == 0L)
      stop("No requested `chroms` overlap obj@loops' chromosomes")
  }
  loop_idx_by_chr <- split(seq_along(loop_chr), loop_chr)

  # ---- worker: normalize one chromosome, return (consensus_rows, values matrix) ----
  per_chr_worker <- function(chr) {
    rows <- loop_idx_by_chr[[chr]]
    want_keys <- loop_key[rows]
    n_samples <- length(samples_vec)

    # 1. Pull whole-genome (i.e. whole-chrom) straw per sample
    sparse_list <- lapply(seq_len(n_samples), function(si) {
      dat <- tryCatch(
        strawr::straw("NONE", hic_paths[si], chr, chr, "BP", resolution),
        error = function(e) NULL
      )
      if (is.null(dat) || nrow(dat) == 0L) {
        return(data.frame(x = integer(0), y = integer(0),
                          counts = numeric(0)))
      }
      # canonicalize x <= y
      sw <- dat$x > dat$y
      if (any(sw)) { t1 <- dat$x[sw]; dat$x[sw] <- dat$y[sw]; dat$y[sw] <- t1 }
      dat
    })

    # 1b. Per-sample cis contribution from this chrom (sum of raw counts).
    # Computed here as a side-effect of the straw pulls we're already doing.
    chrom_cis <- vapply(sparse_list,
                        function(d) sum(as.numeric(d$counts)),
                        numeric(1L))

    # 2. Method-specific normalization
    if (method == "multiHiCcompare") {
      norm_vals <- .normalize_multiHiCcompare(sparse_list, chr, span)
    } else {
      # Filter BNBC emission to consensus (binA, binB) pairs on this chrom.
      # QN math is still computed over the full band; only the final triplet
      # emission is restricted, which keeps the dcast wide-pivot from
      # blowing past R's integer.max on large chroms with many samples.
      keep_pairs_chr <- cbind(loop_binA[rows], loop_binB[rows])
      cis_arg <- if (isTRUE(depth_correct)) obj@cis else NULL
      norm_vals <- .normalize_BNBC(sparse_list, chr,
                                   keep_pairs = keep_pairs_chr,
                                   cis = cis_arg)
    }
    # norm_vals: data.table(chr, region1, region2, IF_1, ..., IF_n)

    # 3. Index consensus loops on this chrom into the normalized table
    norm_keys <- paste(norm_vals$chr, norm_vals$region1, norm_vals$region2, sep = ":")
    m <- match(want_keys, norm_keys)
    if_cols <- grep("^IF_", names(norm_vals), value = TRUE)
    chrom_mat <- matrix(0, nrow = length(rows), ncol = n_samples,
                        dimnames = list(NULL, samples_vec))
    for (j in seq_len(n_samples)) {
      v <- norm_vals[[if_cols[j]]]
      vals <- ifelse(is.na(m), 0, v[m])
      vals[is.na(vals)] <- 0
      chrom_mat[, j] <- vals
    }
    list(rows = rows, mat = chrom_mat, cis = chrom_cis)
  }

  chrom_results <- BiocParallel::bplapply(chroms_to_run, per_chr_worker, BPPARAM = BPPARAM)

  # ---- Assemble per-chrom slices into the full (loops x samples) phenotype ----
  # Rows for chroms we did NOT process are left as NA, so a partial run produces
  # an obj that can be combined with sibling per-chrom runs via
  # loopQTL_combine_normalized().
  pheno <- matrix(NA_real_, nrow = length(loops_gi), ncol = length(samples_vec),
                  dimnames = list(NULL, samples_vec))
  processed_rows <- integer(0)
  cis_accum <- setNames(numeric(length(samples_vec)), samples_vec)
  for (r in chrom_results) {
    pheno[r$rows, ] <- r$mat
    processed_rows <- c(processed_rows, r$rows)
    cis_accum <- cis_accum + r$cis
  }

  # Record cis only if the user did not supply it at setup.
  # If user provided cis (obj@meta$cis_source == "user_provided"), leave it
  # alone — their value is genome-wide and authoritative.
  user_provided_cis <- isTRUE(identical(obj@meta$cis_source, "user_provided"))
  if (!user_provided_cis) {
    obj@cis <- cis_accum
    obj@meta$cis_source     <- "computed_from_hic"
    obj@meta$cis_chroms_in  <- chroms_to_run
  }

  # No RINT here — that's run_QTL()'s job (so obj@phenotype always stores the
  # pre-RINT BNBC/fastlo output, available for QC/comparison).
  obj@phenotype <- pheno
  obj@meta$normalize_method   <- method
  obj@meta$normalize_order    <- method
  obj@meta$chroms_normalized  <- chroms_to_run
  methods::validObject(obj)
  obj
}

# -------------------------------------------------------------------------- #
# Internal: multiHiCcompare cyclic loess (fast) for one chrom across samples
# -------------------------------------------------------------------------- #
.normalize_multiHiCcompare <- function(sparse_list, chr, span) {
  n_samples <- length(sparse_list)
  # multiHiCcompare expects data.frames with columns: chr (int or chr), region1, region2, IF
  dl <- lapply(sparse_list, function(d) {
    if (nrow(d) == 0L) {
      return(data.frame(chr = integer(0),
                        region1 = integer(0),
                        region2 = integer(0),
                        IF = numeric(0)))
    }
    data.frame(
      chr     = .chr_to_int(chr),
      region1 = as.integer(d$x),
      region2 = as.integer(d$y),
      IF      = as.numeric(d$counts)
    )
  })

  # If every sample is empty for this chrom, return an empty merged table
  if (all(vapply(dl, nrow, integer(1L)) == 0L)) {
    out <- data.table::data.table(
      chr = character(0), region1 = integer(0), region2 = integer(0)
    )
    for (j in seq_len(n_samples))
      out[[paste0("IF_", j)]] <- numeric(0)
    return(out)
  }

  hicexp <- multiHiCcompare::make_hicexp(
    data_list = dl,
    groups    = rep(1L, n_samples),
    remove.regions = NULL,
    zero.p    = 1,        # don't drop bin-pairs on % zeros — we want raw coverage
    A.min     = 0
  )
  # fastlo is multiHiCcompare's fast loess normalization (distinct from cyclic_loess).
  # `parallel = FALSE` here because we parallelize one level up, by chromosome.
  hicexp <- multiHiCcompare::fastlo(
    hicexp, span = span, parallel = FALSE, verbose = FALSE
  )
  tab <- multiHiCcompare::hic_table(hicexp)
  # `tab` has columns: chr, region1, region2, D, IF1, IF2, ...
  # Rename IFj -> IF_j and ensure chr is character matching `chr` we expect
  if_cols_in <- grep("^IF\\d+$", names(tab), value = TRUE)
  data.table::setnames(tab, if_cols_in,
                       paste0("IF_", seq_along(if_cols_in)))
  tab$chr <- as.character(chr)
  tab[, c("chr", "region1", "region2", paste0("IF_", seq_along(if_cols_in))), with = FALSE]
}

# -------------------------------------------------------------------------- #
# Internal: BNBC band-wise quantile normalization, one chrom across samples.
# Skips ComBat. Sparse triplet implementation (no dense per-sample matrices).
# QN math is delegated to preprocessCore::normalize.quantiles(), the canonical
# implementation used by BNBC and limma.
#
# Memory profile (chr1 / 50-100 samples / 10 kb): ~25-50 GB peak.
#
# Algorithm:
#   1. For each sample, store straw output as a triplet (i, j, x, d=j-i)
#      data.table in upper-triangle convention. Skip zeros entirely.
#   2. For each genomic distance band d (in bp), find the union of bin-pair
#      positions where any sample has a non-zero contact. Build a small
#      dense L_d x n_samples block with 0-fill for missing samples.
#   3. Quantile-normalize that block across samples via preprocessCore.
#   4. Emit non-zero QN'd cells as new triplets, OPTIONALLY filtered to a
#      user-supplied set of (i, j) keep_pairs. QN math is always over the
#      full band; the filter only restricts what we materialize for the
#      downstream wide pivot.
#   5. Pivot the combined triplets to wide form (one IF column per sample).
#
# Why keep_pairs matters: without it, the long->wide pivot (dcast) tries to
# build a grid of all genome-wide non-zero (i, j) x all samples, which can
# exceed R's integer.max for large chroms with many samples (e.g. chr3 /
# 50 samples / 10 kb gives ~6.5e9 cells in the CJ). Restricting emission
# to the consensus-loop (i, j) pairs that the caller will use shrinks this
# to a few thousand rows per chrom.
#
# Correctness: bin-pairs that are zero in every sample are absent from all
# triplets and stay absent (their QN'd value would be the row-1 reference,
# which equals mean(per-band minima); since most band positions are all-zero
# for chr-scale Hi-C, that reference is ~0). Bin-pairs non-zero in at least
# one sample are processed in full.
# -------------------------------------------------------------------------- #
.normalize_BNBC <- function(sparse_list, chr, keep_pairs = NULL, cis = NULL) {
  n_samples <- length(sparse_list)
  empty_out <- function() {
    out <- data.table::data.table(
      chr = character(0), region1 = integer(0), region2 = integer(0))
    for (k in seq_len(n_samples)) out[[paste0("IF_", k)]] <- numeric(0)
    out
  }
  if (n_samples == 0L) return(empty_out())

  # Pre-compute a hash set of keep_pairs for O(1) lookup per-band.
  # keep_pairs is a 2-column matrix/data.frame of (i, j) in bp.
  keep_set <- NULL
  if (!is.null(keep_pairs)) {
    if (NROW(keep_pairs) == 0L) return(empty_out())
    kp_i <- as.integer(keep_pairs[, 1L])
    kp_j <- as.integer(keep_pairs[, 2L])
    swap <- kp_i > kp_j
    if (any(swap)) { t1 <- kp_i[swap]; kp_i[swap] <- kp_j[swap]; kp_j[swap] <- t1 }
    keep_set <- paste(kp_i, kp_j, sep = "_")
  }

  # ---- 0. logCPM depth correction (exact bnbc formula) ----
  # bnbc/R/contactGroup_utils.R::logCPM:
  #   libs <- (librarySize + 1) / 1e6
  #   c.   <- abs(1/min(libs) - 1/max(libs))      # stability offset
  #   log((count + 0.5) / libs + c.)
  # natural log, not log2; (count + 0.5) not (count + 1).
  # When cis is NULL we fall back to log2(count + 1) (no depth correction).
  if (!is.null(cis)) {
    libs <- (as.numeric(cis) + 1) / 1e6
    if (any(!is.finite(libs)) || any(libs <= 0))
      stop("Invalid cis values: must be positive and finite")
    c_offset <- abs(1/min(libs) - 1/max(libs))
    transform_sample <- function(counts, s) {
      log((as.numeric(counts) + 0.5) / libs[s] + c_offset)
    }
  } else {
    transform_sample <- function(counts, s) {
      log2(as.numeric(counts) + 1)
    }
  }

  # ---- 1. per-sample triplets, depth-corrected log, upper-triangle, keyed by band ----
  per_sample_dt <- vector("list", n_samples)
  for (s in seq_len(n_samples)) {
    d <- sparse_list[[s]]
    if (is.null(d) || nrow(d) == 0L) {
      per_sample_dt[[s]] <- data.table::data.table(
        i = integer(0), j = integer(0),
        x = numeric(0), d = integer(0))
      data.table::setkey(per_sample_dt[[s]], d, i)
      next
    }
    ii <- as.integer(d$x)
    jj <- as.integer(d$y)
    swap <- ii > jj
    if (any(swap)) { t1 <- ii[swap]; ii[swap] <- jj[swap]; jj[swap] <- t1 }
    dt <- data.table::data.table(
      i = ii, j = jj,
      x = transform_sample(d$counts, s),
      d = jj - ii
    )
    data.table::setkey(dt, d, i)
    per_sample_dt[[s]] <- dt
  }

  # ---- 2. set of distinct bands across all samples ----
  all_bands <- sort(unique(unlist(
    lapply(per_sample_dt, function(dt) unique(dt$d))
  )))
  if (length(all_bands) == 0L) return(empty_out())

  # ---- 3. per-band quantile normalization ----
  output_chunks <- vector("list", length(all_bands))
  for (k in seq_along(all_bands)) {
    b <- all_bands[k]

    # Per-sample (i, x) for this band via keyed lookup
    samp_data <- lapply(per_sample_dt, function(dt) {
      dt[.(b), .(i, x), on = "d", nomatch = 0L]
    })

    # Union of active i positions for this band
    active_i <- sort(unique(unlist(lapply(samp_data, function(sd) sd$i))))
    L <- length(active_i)
    if (L == 0L) next
    # Note: L = 1 is OK — preprocessCore returns all samples equal to the mean.
    # That degenerate cell will have zero variance downstream and be skipped.

    # Build L x n_samples dense block. The "missing" fill is the per-sample
    # logCPM transform of raw-count=0 (not literally 0) — because not having
    # a triplet at (i, i+b) means the underlying raw count was 0, and
    # log((0+0.5)/libs[s] + c_offset) is a small negative value, not 0.
    zero_fill <- vapply(seq_len(n_samples), function(s) transform_sample(0, s),
                        numeric(1L))
    block <- matrix(zero_fill, nrow = L, ncol = n_samples, byrow = TRUE)
    for (s in seq_len(n_samples)) {
      sd <- samp_data[[s]]
      if (nrow(sd) > 0L) {
        row_idx <- match(sd$i, active_i)
        block[row_idx, s] <- sd$x
      }
    }

    qn <- .quantile_normalize(block)

    # Filter QN output rows to keep_pairs (if requested) BEFORE materializing
    # the triplet emission. This is the fix that prevents the dcast CJ
    # overflow on large chroms with many samples.
    if (!is.null(keep_set)) {
      band_keys <- paste(active_i, active_i + b, sep = "_")
      keep_rows <- which(band_keys %in% keep_set)
      if (length(keep_rows) == 0L) next
      qn <- qn[keep_rows, , drop = FALSE]
      active_i_em <- active_i[keep_rows]
    } else {
      active_i_em <- active_i
    }

    # Emit QN'd cells as triplets. With logCPM, "zero contact" maps to a
    # per-sample baseline value (small but informative), so we emit ALL
    # cells in the keep set. With log2(x+1) (no depth correction), genuine
    # zeros stay zero and we skip them to save downstream memory.
    if (is.null(cis)) {
      nz_mask <- qn != 0
      if (!any(nz_mask)) next
      idxs <- which(nz_mask, arr.ind = TRUE)
      vals <- qn[nz_mask]
    } else {
      # Emit every cell (all are informative under logCPM)
      idxs <- which(matrix(TRUE, nrow = nrow(qn), ncol = ncol(qn)), arr.ind = TRUE)
      vals <- as.numeric(qn)
    }
    output_chunks[[k]] <- data.table::data.table(
      s = idxs[, 2L],
      i = active_i_em[idxs[, 1L]],
      j = active_i_em[idxs[, 1L]] + b,
      x = vals
    )
  }

  # Free per-sample input triplets before assembling the wide output
  rm(per_sample_dt)
  invisible(gc(verbose = FALSE))

  # ---- 4. combine all output chunks into long form ----
  out_long <- data.table::rbindlist(output_chunks, fill = TRUE)
  rm(output_chunks)
  if (nrow(out_long) == 0L) return(empty_out())

  # ---- 5. pivot wide on (i, j) -> one IF column per sample ----
  wide <- data.table::dcast(out_long, i + j ~ s, value.var = "x", fill = 0)
  rm(out_long); invisible(gc(verbose = FALSE))
  data.table::setnames(wide, c("i", "j"), c("region1", "region2"))

  # Ensure every sample has an IF column (some may be all-absent on this chrom)
  for (s in seq_len(n_samples)) {
    nm_in  <- as.character(s)
    nm_out <- paste0("IF_", s)
    if (nm_in %in% names(wide)) {
      data.table::setnames(wide, nm_in, nm_out)
    } else {
      wide[, (nm_out) := 0]
    }
  }

  wide[, chr := as.character(chr)]
  wide[, c("chr", "region1", "region2",
           paste0("IF_", seq_len(n_samples))), with = FALSE]
}

# Quantile normalization across columns (samples) of an L x N matrix.
# Thin wrapper around preprocessCore::normalize.quantiles(), which is the
# canonical Bioconductor QN implementation used by BNBC, limma, et al.
# (Bolstad et al. 2003). Single-column input is returned unchanged.
.quantile_normalize <- function(mat) {
  if (ncol(mat) < 2L) return(mat)
  out <- preprocessCore::normalize.quantiles(mat, copy = TRUE)
  dimnames(out) <- dimnames(mat)
  out
}

# Convert chrom name to int form expected by multiHiCcompare ("chr1" -> 1, "X" -> 23, etc.)
.chr_to_int <- function(chr) {
  ch <- sub("^chr", "", as.character(chr))
  v  <- suppressWarnings(as.integer(ch))
  if (!is.na(v)) return(v)
  switch(ch, "X" = 23L, "Y" = 24L, "M" = 25L, "MT" = 25L,
         # fallback: hash to a stable positive int
         abs(as.integer(charToRaw(ch)[1L])) + 100L)
}
