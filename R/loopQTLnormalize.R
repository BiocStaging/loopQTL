#' Normalize loopQTL contacts and write the per-loop phenotype matrix
#'
#' For each chromosome, pulls whole-chrom `(x, y, counts)` from every .hic
#' via `strawr::straw("NONE", ...)`, applies per-sample logCPM depth
#' correction (`log((count + 0.5) / libs + c_offset)`, requires `obj@cis`),
#' then band-wise quantile-normalizes each off-diagonal band across samples
#' via `preprocessCore::normalize.quantiles()` (Bolstad et al. 2003),
#' following Fletez-Brant et al. 2021. Normalized values are looked up at
#' consensus bin-pairs. Sparse triplet implementation; the `bnbc`
#' Bioconductor package is not used. Per-chromosome work is dispatched via
#' `BiocParallel`.
#'
#' **No RINT is applied here** - the phenotype slot stores normalize's output
#' directly. RINT is applied just before regression in `runQTL()`
#' (`normal = TRUE` default).
#'
#' @param obj A [loopQTL-class] from [loopQTLsetup()].
#' @param BPPARAM A [BiocParallel::BiocParallelParam] for parallelization
#'   over chromosomes. Default [BiocParallel::SerialParam()].
#' @param save_windows Logical. If `TRUE`, also save a local contact submatrix
#'   around each in-loop SNP for downstream per-pixel SNP-correlation maps
#'   (Shi et al. 2025). Values are sliced from the band-QN'd per-chrom matrix
#'   in the same pass that fills the phenotype, so they inherit the same
#'   depth correction + cross-sample QN. Purely descriptive output; does not
#'   feed back into `obj@phenotype` or `runQTL()`.
#' @param windows_out Required when `save_windows = TRUE`. Path to a `.rds`
#'   file the per-SNP windows are written to (xz-compressed). The file
#'   contains `list(cells, cell_coords, snp_index, samples, resolution)`;
#'   see Details.
#' @param window Minimum window extent in bp around each in-loop SNP.
#'   Default 1e6. The actual window is
#'   `max(window, span(SNP_pos U all partner anchors of loops the SNP
#'   participates in))`, centered on the bounding interval. Guarantees the
#'   partner anchor stays in-frame even for ~2 Mb Mustache loops. Ignored
#'   when `save_windows = FALSE`.
#' @param snp_window Bp flank added to each consensus anchor when selecting
#'   in-loop SNPs from the VCF. Default 0 (anchor-bin SNPs only - matches
#'   `runQTL()`'s default). Ignored when `save_windows = FALSE`.
#'
#' @return The input object with `phenotype` filled. When `save_windows = TRUE`,
#'   `meta$windows_out` records the ABSOLUTE file path (via `normalizePath()`,
#'   so downstream tools like `loopQTLplotSnpWindow()` can find the file
#'   regardless of the working directory at read time), `meta$n_windows_saved`
#'   the SNP count.
#'
#' @details
#' **Saved windows file format (`.rds`):** a list with
#' \describe{
#'   \item{`cells`}{numeric matrix `samples x N_cells`. One column per
#'     `(chrom, bin_i, bin_j)` cell in the union of every in-loop SNP's
#'     window box. Values are post-band-QN, pre-RINT (contact-scale).
#'     `NA -> 0` fill matches loop extraction. Sample axis follows `samples`.}
#'   \item{`cells_raw`}{numeric matrix `samples x N_cells`, same shape and
#'     column order as `cells`, but holding the RAW straw counts (before
#'     any depth correction or QN). Divide row-wise by `obj@cis / 1e6`
#'     to get per-sample CPM, or compute your own transform.}
#'   \item{`cell_coords`}{data.frame with `chrom`, `bin_i`, `bin_j` per
#'     column of `cells` / `cells_raw`.}
#'   \item{`snp_index`}{data.frame with one row per SNP: `snp_id`, `chrom`,
#'     `pos`, `ref`, `alt`, `bin_lo`, `bin_hi`, `n_loops_anchored`, plus a
#'     list-column `col_indices` of integer column indices into `cells`.}
#'   \item{`samples`}{character vector matching the rows of `cells`.}
#'   \item{`resolution`}{integer; bin size in bp.}
#' }
#'
#' @export
loopQTLnormalize <- function(obj,
                              BPPARAM       = BiocParallel::SerialParam(),
                              save_windows  = FALSE,
                              windows_out   = NULL,
                              window        = 1e6,
                              snp_window    = 0L) {
  stopifnot(methods::is(obj, "loopQTL"))

  # save_windows validation: path required, VCF must be set
  if (isTRUE(save_windows)) {
    if (is.null(windows_out) || !is.character(windows_out) ||
        length(windows_out) != 1L || !nzchar(windows_out))
      stop("save_windows = TRUE requires `windows_out`: a single .rds path")
    stopifnot(is.numeric(window),     length(window)     == 1L, window     > 0,
              is.numeric(snp_window), length(snp_window) == 1L, snp_window >= 0)
    window     <- as.numeric(window)
    snp_window <- as.integer(snp_window)
    # Best-effort writable check (touch the file)
    ok_touch <- tryCatch({
      con <- file(windows_out, open = "wb"); close(con); file.remove(windows_out); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok_touch)
      stop("Cannot write to windows_out path: ", windows_out)
    if (is.null(obj@vcf_path) || is.na(obj@vcf_path) || !nzchar(obj@vcf_path))
      stop("obj@vcf_path is not set; cannot read SNPs for windows")
  } else if (!is.null(windows_out)) {
    warning("`windows_out` is ignored when save_windows = FALSE")
  }

  # BNBC's logCPM depth correction requires cis
  if (length(obj@cis) == 0L) {
    stop("Normalize requires per-sample cis depth in obj@cis, but it is ",
         "empty. Options:\n",
         "  1. Include a `cis` column in the samplesheet at loopQTLsetup(), or\n",
         "  2. Re-run loopQTLsetup() (cis is auto-computed from .hic when no\n",
         "     `cis` column is provided), or\n",
         "  3. Manually attach: obj@cis <- setNames(<numeric>, obj@samples)")
  }
  if (length(obj@cis) != length(obj@samples) ||
      !identical(names(obj@cis), obj@samples)) {
    stop("obj@cis must be a named numeric vector with one entry per sample, ",
         "names matching obj@samples")
  }

  samples_vec <- obj@samples
  resolution  <- obj@resolution
  loops_gi    <- obj@loops

  # Per-sample per-chrom hic paths from meta (we re-derive from setup invariants)
  if (is.null(obj@meta$hic_paths)) {
    stop("obj@meta$hic_paths missing; loopQTLsetup must record the per-sample .hic paths")
  }
  hic_paths <- obj@meta$hic_paths

  # Consensus loops keyed by chrom + canonical (binA <= binB)
  anchorA <- InteractionSet::anchors(loops_gi, type = "first")
  anchorB <- InteractionSet::anchors(loops_gi, type = "second")
  loop_chr  <- as.character(GenomicRanges::seqnames(anchorA))
  # Floor each anchor's midpoint to the bin containing it - straw returns
  # bin-based contacts, so we look up at bin granularity even when the consensus
  # key is exact-coord.
  loop_binA <- .bin_anchor(GenomicRanges::start(anchorA),
                           GenomicRanges::end(anchorA), resolution)
  loop_binB <- .bin_anchor(GenomicRanges::start(anchorB),
                           GenomicRanges::end(anchorB), resolution)
  swap <- loop_binA > loop_binB
  if (any(swap)) {
    tmp <- loop_binA[swap]; loop_binA[swap] <- loop_binB[swap]; loop_binB[swap] <- tmp
  }
  loop_key <- paste(loop_chr, loop_binA, loop_binB, sep = ":")

  chroms_to_run <- unique(loop_chr)
  loop_idx_by_chr <- split(seq_along(loop_chr), loop_chr)

  # Per-sample chrom-naming normalizers so BNBC straw calls translate the
  # anchor-side chrom ("chr10") to whatever each .hic file uses ("10" or
  # "chr10"). Without this, BEDPE-vs-hic naming mismatch makes straw
  # silently return empty data -> all-zero phenotype -> all-NA after RINT.
  chrom_norm_BNBC <- lapply(hic_paths, function(p) {
    .make_chrom_normalizer(strawr::readHicChroms(p)$name)
  })

  # save_windows: precompute SNP -> (chrom, cells to slice)
  windows_ctx <- NULL
  if (isTRUE(save_windows)) {
    windows_ctx <- .compute_snp_windows(
      vcf_path   = obj@vcf_path,
      anchorA    = anchorA,
      anchorB    = anchorB,
      snp_window = snp_window,
      window     = window,
      resolution = resolution
    )
    if (is.null(windows_ctx)) {
      message("[", format(Sys.time()), "] save_windows: nothing to save; ",
              "continuing without window output")
      save_windows <- FALSE
    } else {
      # Sanity check: every SNP's chrom must be in our chroms_to_run
      win_chroms <- names(windows_ctx$cells_by_chrom)
      missing_chr <- setdiff(win_chroms, chroms_to_run)
      if (length(missing_chr) > 0L)
        stop("save_windows: SNP chrom(s) not in normalize's chroms_to_run: ",
             paste(missing_chr, collapse = ", "))
      n_total_cells <- sum(vapply(windows_ctx$cells_by_chrom,
                                  nrow, integer(1L)))
      est_mb <- n_total_cells * length(samples_vec) * 8 / 1024^2
      message(sprintf(
        "[%s] save_windows: %d SNP(s), %d unique cells across %d chrom(s); estimated ~%.1f MB in memory (8 bytes x cells x samples)",
        format(Sys.time()), nrow(windows_ctx$snps_df),
        n_total_cells, length(win_chroms), est_mb))
    }
  }

  # worker: normalize one chromosome, return (consensus_rows, values matrix)
  per_chr_worker <- function(chr) {
    rows <- loop_idx_by_chr[[chr]]
    want_keys <- loop_key[rows]
    n_samples <- length(samples_vec)

    # 1. Pull whole-genome (i.e. whole-chrom) straw per sample
    sparse_list <- lapply(seq_len(n_samples), function(si) {
      chr_in_hic <- chrom_norm_BNBC[[si]](chr)
      dat <- tryCatch(
        strawr::straw("NONE", hic_paths[si], chr_in_hic, chr_in_hic,
                      "BP", resolution),
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

    # 2. BNBC normalization. Filter emission to consensus (binA, binB) pairs
    # on this chrom; QN math is still computed over the full band, only the
    # final triplet emission is restricted, which keeps the dcast wide-pivot
    # from blowing past R's integer.max on large chroms with many samples.
    # When save_windows is on, also include the SNP window cells in the keep
    # set so a single BNBC pass emits both consensus + window cells.
    keep_pairs_chr <- cbind(loop_binA[rows], loop_binB[rows])
    window_cells_chr <- NULL
    if (isTRUE(save_windows) && !is.null(windows_ctx) &&
        !is.null(windows_ctx$cells_by_chrom[[chr]])) {
      window_cells_chr <- windows_ctx$cells_by_chrom[[chr]]
      combined <- rbind(
        keep_pairs_chr,
        cbind(window_cells_chr$bin_i, window_cells_chr$bin_j)
      )
      keep_pairs_pass <- unique(combined)
    } else {
      keep_pairs_pass <- keep_pairs_chr
    }
    cis_arg <- obj@cis
    norm_vals <- .normalize_BNBC(sparse_list, chr,
                                 keep_pairs = keep_pairs_pass,
                                 cis = cis_arg)
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

    # 4. (save_windows) Index window cells into two (samples x n_cells) matrices:
    #    - windows_mat:     post-BNBC normalized values (from norm_vals)
    #    - windows_raw_mat: raw straw counts at same cells (from sparse_list)
    windows_mat     <- NULL
    windows_raw_mat <- NULL
    if (!is.null(window_cells_chr)) {
      win_keys <- paste(chr, window_cells_chr$bin_i, window_cells_chr$bin_j,
                        sep = ":")
      mw <- match(win_keys, norm_keys)
      windows_mat <- matrix(0, nrow = n_samples, ncol = nrow(window_cells_chr),
                            dimnames = list(samples_vec, NULL))
      for (j in seq_len(n_samples)) {
        v <- norm_vals[[if_cols[j]]]
        vals <- ifelse(is.na(mw), 0, v[mw])
        vals[is.na(vals)] <- 0
        windows_mat[j, ] <- vals
      }

      # Raw counts: match window cells against each sample's straw sparse output
      win_key_local <- paste(window_cells_chr$bin_i, window_cells_chr$bin_j,
                             sep = ":")
      windows_raw_mat <- matrix(0, nrow = n_samples,
                                ncol = nrow(window_cells_chr),
                                dimnames = list(samples_vec, NULL))
      for (s in seq_len(n_samples)) {
        dat <- sparse_list[[s]]
        if (is.null(dat) || nrow(dat) == 0L) next
        dat_key <- paste(dat$x, dat$y, sep = ":")
        m <- match(win_key_local, dat_key)
        windows_raw_mat[s, ] <- ifelse(is.na(m), 0, dat$counts[m])
      }
    }
    list(rows = rows, mat = chrom_mat,
         windows_mat     = windows_mat,
         windows_raw_mat = windows_raw_mat,
         window_cells_chr = window_cells_chr)
  }

  chrom_results <- BiocParallel::bplapply(chroms_to_run, per_chr_worker, BPPARAM = BPPARAM)
  names(chrom_results) <- chroms_to_run

  # Assemble per-chrom slices into the full (loops x samples) phenotype
  pheno <- matrix(NA_real_, nrow = length(loops_gi), ncol = length(samples_vec),
                  dimnames = list(NULL, samples_vec))
  for (r in chrom_results) {
    pheno[r$rows, ] <- r$mat
  }

  # cis is set authoritatively at setup time and is genome-wide; normalize
  # does not touch it.

  # No RINT here - that's runQTL()'s job (so obj@phenotype always stores the
  # pre-RINT BNBC output, available for QC/comparison).
  obj@phenotype <- pheno

  # Stitch per-chrom window outputs and write .rds
  if (isTRUE(save_windows) && !is.null(windows_ctx)) {
    # Per-chrom column offsets into the genome-wide cells matrix
    chr_order <- intersect(chroms_to_run, names(windows_ctx$cells_by_chrom))
    n_cells_per_chr <- vapply(chr_order, function(chr) {
      r <- chrom_results[[chr]]
      if (is.null(r$windows_mat)) 0L else ncol(r$windows_mat)
    }, integer(1L))
    if (sum(n_cells_per_chr) == 0L) {
      message("[", format(Sys.time()), "] save_windows: every per-chrom ",
              "window slice came back empty; nothing to write")
    } else {
      total_cells <- sum(n_cells_per_chr)
      cells     <- matrix(0, nrow = length(samples_vec), ncol = total_cells,
                          dimnames = list(samples_vec, NULL))
      cells_raw <- matrix(0, nrow = length(samples_vec), ncol = total_cells,
                          dimnames = list(samples_vec, NULL))
      cell_coords <- data.frame(
        chrom = character(total_cells),
        bin_i = integer(total_cells),
        bin_j = integer(total_cells),
        stringsAsFactors = FALSE
      )
      offset <- 0L
      chrom_col_offset <- setNames(integer(length(chr_order)), chr_order)
      for (chr in chr_order) {
        r <- chrom_results[[chr]]
        if (is.null(r$windows_mat) || ncol(r$windows_mat) == 0L) {
          chrom_col_offset[chr] <- NA_integer_
          next
        }
        ncol_chr <- ncol(r$windows_mat)
        rng <- seq.int(offset + 1L, offset + ncol_chr)
        cells[, rng]     <- r$windows_mat
        cells_raw[, rng] <- r$windows_raw_mat
        cell_coords$chrom[rng] <- chr
        cell_coords$bin_i[rng] <- r$window_cells_chr$bin_i
        cell_coords$bin_j[rng] <- r$window_cells_chr$bin_j
        chrom_col_offset[chr]  <- offset
        offset <- offset + ncol_chr
      }

      # Per-SNP col_indices: local within-chrom indices + that chrom's offset
      snp_chrom <- windows_ctx$snps_df$chrom_orig
      n_snp <- nrow(windows_ctx$snps_df)
      col_indices <- vector("list", n_snp)
      for (i in seq_len(n_snp)) {
        chr <- snp_chrom[i]
        off <- chrom_col_offset[chr]
        loc <- windows_ctx$snp_cells_local[[i]]
        col_indices[[i]] <- if (is.na(off) || length(loc) == 0L) integer(0)
                            else as.integer(loc + off)
      }

      snp_index <- windows_ctx$snps_df
      snp_index$chrom <- snp_index$chrom_orig
      snp_index$chrom_orig <- NULL
      snp_index$chrom_vcf  <- NULL
      snp_index$col_indices <- I(col_indices)
      snp_index <- snp_index[, c("snp_id", "chrom", "pos", "ref", "alt",
                                  "bin_lo", "bin_hi", "n_loops_anchored",
                                  "col_indices")]

      saveRDS(
        list(
          cells       = cells,           # post-BNBC normalized values
          cells_raw   = cells_raw,       # raw straw counts (integer)
          cell_coords = cell_coords,
          snp_index   = snp_index,
          samples     = samples_vec,
          resolution  = resolution
        ),
        file = windows_out,
        compress = "xz"
      )
      message(sprintf(
        "[%s] save_windows: wrote %s (%.1f MB)",
        format(Sys.time()), windows_out,
        file.info(windows_out)$size / 1024^2))

      # Store the ABSOLUTE path so downstream loaders (loopQTLplotSnpWindow,
      # cross-session .rda restores) find the file regardless of the working
      # directory in effect at read time. Falls back to the input if the file
      # somehow can't be resolved (should never happen since we just wrote it).
      obj@meta$windows_out      <- tryCatch(
        normalizePath(windows_out, winslash = "/", mustWork = TRUE),
        error = function(e) windows_out
      )
      obj@meta$n_windows_saved  <- nrow(snp_index)
      obj@meta$n_window_cells   <- total_cells
      obj@meta$window_bp        <- as.numeric(window)
      obj@meta$snp_window_bp    <- as.integer(snp_window)
    }
  }

  methods::validObject(obj)
  obj
}

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

  # 0. logCPM depth correction (exact bnbc formula)
  #    bnbc/R/contactGroup_utils.R::logCPM:
  #      libs <- (librarySize + 1) / 1e6
  #      c.   <- abs(1/min(libs) - 1/max(libs))      # stability offset
  #      log((count + 0.5) / libs + c.)
  #    natural log, not log2; (count + 0.5) not (count + 1).
  #    When cis is NULL we fall back to log2(count + 1) (no depth correction).
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

  # 1. per-sample triplets, depth-corrected log, upper-triangle, keyed by band
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

  # 2. set of distinct bands across all samples
  all_bands <- sort(unique(unlist(
    lapply(per_sample_dt, function(dt) unique(dt$d))
  )))
  if (length(all_bands) == 0L) return(empty_out())

  # 3. per-band quantile normalization
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
    # Note: L = 1 is OK - preprocessCore returns all samples equal to the mean.
    # That degenerate cell will have zero variance downstream and be skipped.

    # Build L x n_samples dense block. The "missing" fill is the per-sample
    # logCPM transform of raw-count=0 (not literally 0) - because not having
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

  # 4. combine all output chunks into long form
  out_long <- data.table::rbindlist(output_chunks, fill = TRUE)
  rm(output_chunks)
  if (nrow(out_long) == 0L) return(empty_out())

  # 5. pivot wide on (i, j) -> one IF column per sample
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

