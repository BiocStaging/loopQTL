#' Normalize loopQTL contacts at whole-genome scale, optional RINT
#'
#' For each chromosome, pulls whole-genome `(x, y, counts)` from every .hic via
#' `strawr::straw("NONE", ...)`, then normalizes across samples:
#'
#'  * `method = "multiHiCcompare"`: builds a per-chrom `hicexp` (groups all = 1)
#'    and runs `multiHiCcompare::fastlo(...)`. `fastlo` is multiHiCcompare's
#'    fast loess normalization, distinct from `cyclic_loess`. We do not use
#'    the differential side of multiHiCcompare.
#'  * `method = "BNBC"`: per-chrom log-transformed contact matrices,
#'    band-wise quantile-normalized across samples. ComBat batch correction is
#'    intentionally skipped.
#'
#' Normalized values are indexed at the consensus loop bin-pairs to build the
#' (loops x samples) phenotype matrix. If `normal = TRUE` (default), a Blom-
#' offset RINT is applied per loop across samples. Zero-variance rows become
#' all-NA phenotype rows and are counted in `meta$n_zero_var`.
#'
#' Per-chromosome work is dispatched via `BiocParallel`. This normalize-time
#' straw pull is independent of (and additional to) the consensus extraction
#' performed at `loopQTL_setup()` time.
#'
#' @param obj A [loopQTL-class] from [loopQTL_setup()].
#' @param method Normalization method. Default `"multiHiCcompare"`.
#' @param normal Apply RINT after normalization. Default `TRUE`.
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
#' @return The input object with `phenotype` filled and `meta$normalize_method`,
#'   `meta$normal`, `meta$n_zero_var` recorded.
#'
#' @export
loopQTL_normalize <- function(obj,
                              method  = c("multiHiCcompare", "BNBC"),
                              normal  = TRUE,
                              span    = 0.7,
                              chroms  = NULL,
                              BPPARAM = BiocParallel::SerialParam()) {
  stopifnot(methods::is(obj, "loopQTL"))
  method <- match.arg(method)

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

    # 2. Method-specific normalization
    if (method == "multiHiCcompare") {
      norm_vals <- .normalize_multiHiCcompare(sparse_list, chr, span)
    } else {
      norm_vals <- .normalize_BNBC(sparse_list, chr)
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
    list(rows = rows, mat = chrom_mat)
  }

  chrom_results <- BiocParallel::bplapply(chroms_to_run, per_chr_worker, BPPARAM = BPPARAM)

  # ---- Assemble per-chrom slices into the full (loops x samples) phenotype ----
  # Rows for chroms we did NOT process are left as NA, so a partial run produces
  # an obj that can be combined with sibling per-chrom runs via
  # loopQTL_combine_normalized().
  pheno <- matrix(NA_real_, nrow = length(loops_gi), ncol = length(samples_vec),
                  dimnames = list(NULL, samples_vec))
  processed_rows <- integer(0)
  for (r in chrom_results) {
    pheno[r$rows, ] <- r$mat
    processed_rows <- c(processed_rows, r$rows)
  }

  # ---- RINT if requested (per-loop; only over the rows we actually filled) ----
  n_zero_var <- 0L
  if (normal && length(processed_rows) > 0L) {
    sub <- pheno[processed_rows, , drop = FALSE]
    row_sd <- apply(sub, 1L, stats::sd)
    zero_var <- !is.finite(row_sd) | row_sd == 0
    n_zero_var <- sum(zero_var)
    rint_sub <- sub
    rint_sub[zero_var, ] <- NA_real_
    if (any(!zero_var)) {
      good <- which(!zero_var)
      rint_sub[good, ] <- t(apply(sub[good, , drop = FALSE], 1L, .rint))
    }
    pheno[processed_rows, ] <- rint_sub
  }

  obj@phenotype <- pheno
  obj@meta$normalize_method   <- method
  obj@meta$normal             <- normal
  obj@meta$n_zero_var         <- n_zero_var
  obj@meta$normalize_order    <- c(method, if (normal) "RINT" else NULL)
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
# Skips ComBat.
# -------------------------------------------------------------------------- #
.normalize_BNBC <- function(sparse_list, chr) {
  n_samples <- length(sparse_list)

  # 1. Determine universe of bins on this chrom across samples
  all_x <- unique(unlist(lapply(sparse_list, function(d) c(d$x, d$y))))
  if (length(all_x) == 0L) {
    out <- data.table::data.table(chr = character(0),
                                  region1 = integer(0),
                                  region2 = integer(0))
    for (j in seq_len(n_samples))
      out[[paste0("IF_", j)]] <- numeric(0)
    return(out)
  }
  all_bins <- sort(unique(all_x))
  bin_to_idx <- setNames(seq_along(all_bins), as.character(all_bins))
  nb <- length(all_bins)

  # 2. Build per-sample dense contact matrices (upper triangle, on this chrom)
  mats <- lapply(sparse_list, function(d) {
    m <- matrix(0, nrow = nb, ncol = nb)
    if (nrow(d) > 0L) {
      i <- bin_to_idx[as.character(d$x)]
      j <- bin_to_idx[as.character(d$y)]
      m[cbind(i, j)] <- d$counts
      m[cbind(j, i)] <- d$counts
    }
    log2(m + 1)
  })

  # 3. Band-wise quantile normalize across samples.
  #    For each off-diagonal band b in 0:(nb-1), collect the b-th superdiagonal
  #    from every sample as a column, run quantile normalization across columns,
  #    write back.
  for (b in 0:(nb - 1L)) {
    idx_row <- seq_len(nb - b)
    idx_col <- idx_row + b
    cells <- cbind(idx_row, idx_col)
    band_mat <- vapply(mats, function(m) m[cells], numeric(length(idx_row)))
    if (length(idx_row) == 1L) {
      # single cell — nothing to qnorm; just keep
      next
    }
    qn <- .quantile_normalize(band_mat)
    for (s in seq_len(n_samples)) {
      mats[[s]][cells] <- qn[, s]
      mats[[s]][cbind(idx_col, idx_row)] <- qn[, s]   # mirror
    }
  }

  # 4. Convert back to long form
  long_list <- vector("list", n_samples)
  for (s in seq_len(n_samples)) {
    m <- mats[[s]]
    idx <- which(upper.tri(m, diag = TRUE) & m != 0, arr.ind = TRUE)
    if (nrow(idx) == 0L) {
      long_list[[s]] <- data.table::data.table(
        region1 = integer(0), region2 = integer(0), val = numeric(0)
      )
      next
    }
    long_list[[s]] <- data.table::data.table(
      region1 = all_bins[idx[, 1L]],
      region2 = all_bins[idx[, 2L]],
      val     = m[idx]
    )
  }

  # 5. Merge by (region1, region2)
  merged <- long_list[[1L]]
  data.table::setnames(merged, "val", "IF_1")
  if (n_samples > 1L) {
    for (s in seq.int(2L, n_samples)) {
      x <- long_list[[s]]
      data.table::setnames(x, "val", paste0("IF_", s))
      merged <- merge(merged, x, by = c("region1", "region2"), all = TRUE)
    }
  }
  for (cc in grep("^IF_", names(merged), value = TRUE)) {
    merged[is.na(merged[[cc]]), (cc) := 0]
  }
  merged$chr <- as.character(chr)
  merged[, c("chr", "region1", "region2",
             paste0("IF_", seq_len(n_samples))), with = FALSE]
}

# Quantile normalization across columns (samples) of a matrix.
# Each column is treated as a sample's vector; each row's mean rank value
# replaces the original.
.quantile_normalize <- function(mat) {
  if (ncol(mat) < 2L) return(mat)
  ranks  <- apply(mat, 2L, rank, ties.method = "average")
  sorted <- apply(mat, 2L, sort)
  row_means <- rowMeans(sorted)
  out <- mat
  for (j in seq_len(ncol(mat))) {
    out[, j] <- row_means[as.integer(ranks[, j])]
  }
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
