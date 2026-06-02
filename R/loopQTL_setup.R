#' Set up a loopQTL analysis
#'
#' Builds the bin-matched consensus loop set, validates per-sample `.hic`
#' metadata and the VCF header, and records the per-sample `.hic` paths
#' needed by [loopQTL_normalize()]. When the samplesheet does not provide a
#' `cis` column, per-sample cis depth is unconditionally auto-computed from
#' the `.hic` files via a one-pass straw aggregation (sums of counts per
#' `(sample, chrom)`, no consensus indexing) — needed by the BNBC logCPM
#' step downstream. If you precompute `cis` per sample yourself (e.g. from
#' a Hi-C QC report) and put it in the samplesheet, that value is used
#' as-is and the .hic pass is skipped.
#'
#' @param samples data.frame with columns `Sample`, `hic`, `loops`. Optional
#'   `cis` column: if present, those values become `obj@cis` directly. If
#'   absent, cis is summed from the `.hic` files here at setup time.
#' @param vcf Path to a bgzipped + tabix-indexed multi-sample VCF.
#' @param resolution Analysis bin size in bp (must exist in every .hic).
#' @param min_samples Keep consensus loops appearing in at least this many
#'   samples (default 2).
#' @param BPPARAM A [BiocParallel::BiocParallelParam] that parallelizes the
#'   cis-computing straw pulls. Default [BiocParallel::SerialParam()].
#'
#' @return A [loopQTL-class] object with `samples`, `loops`, `vcf_path`,
#'   `resolution`, `cis` (if provided or computed), and `meta` populated.
#'   `contacts` and `phenotype` are empty placeholders until normalize fills
#'   `phenotype`.
#'
#' @export
loopQTL_setup <- function(samples,
                          vcf,
                          resolution,
                          min_samples = 2L,
                          BPPARAM     = BiocParallel::SerialParam()) {

  stopifnot(
    is.data.frame(samples),
    all(c("Sample", "hic", "loops") %in% names(samples)),
    is.character(vcf), length(vcf) == 1L,
    is.numeric(resolution), length(resolution) == 1L, resolution > 0,
    is.numeric(min_samples), length(min_samples) == 1L, min_samples >= 1
  )
  resolution  <- as.integer(resolution)
  min_samples <- as.integer(min_samples)

  sample_ids <- as.character(samples$Sample)
  if (anyDuplicated(sample_ids))
    stop("Sample IDs must be unique")
  hic_paths  <- as.character(samples$hic)
  loop_paths <- as.character(samples$loops)
  user_cis   <- if ("cis" %in% names(samples)) as.numeric(samples$cis) else NULL

  for (p in c(hic_paths, loop_paths))
    if (!file.exists(p)) stop("File not found: ", p)

  # ---- VCF validation (header only) ----
  if (!file.exists(vcf)) stop("VCF not found: ", vcf)
  tbi <- paste0(vcf, ".tbi")
  if (!file.exists(tbi))
    stop("Tabix index not found: ", tbi,
         "\nPlease bgzip + tabix-index your VCF (e.g. `tabix -p vcf my.vcf.gz`).")

  vcf_samples <- VariantAnnotation::samples(
    VariantAnnotation::scanVcfHeader(vcf)
  )
  missing_in_vcf <- setdiff(sample_ids, vcf_samples)
  if (length(missing_in_vcf) > 0L)
    stop("Samples in `samples` not present in VCF: ",
         paste(missing_in_vcf, collapse = ", "))

  # ---- Per-.hic metadata + resolution / NONE check (header reads only) ----
  # Also build per-sample chromosome-naming normalizers so we can map a
  # user-style chrom string (e.g. "1") to whatever this .hic file uses
  # (e.g. "chr1") when summing cis.
  chrom_norm_per_sample <- vector("list", length(hic_paths))
  hic_chroms_per_sample <- vector("list", length(hic_paths))
  for (i in seq_along(hic_paths)) {
    chroms_i <- strawr::readHicChroms(hic_paths[i])$name
    bp_res   <- strawr::readHicBpResolutions(hic_paths[i])
    norms    <- strawr::readHicNormTypes(hic_paths[i])
    if (!(resolution %in% bp_res))
      stop("Resolution ", resolution, " not available in ", hic_paths[i],
           " (available: ", paste(bp_res, collapse = ", "), ")")
    if (!("NONE" %in% norms))
      stop("Normalization 'NONE' not available in ", hic_paths[i])
    chrom_norm_per_sample[[i]] <- .make_chrom_normalizer(chroms_i)
    hic_chroms_per_sample[[i]] <- chroms_i
  }

  # ---- Read and partition loops ----
  per_sample_loops <- lapply(loop_paths, .read_loops_bedpe)
  n_trans_per_sample <- vapply(per_sample_loops,
                               function(d) sum(d$chrA != d$chrB),
                               integer(1L))
  per_sample_loops <- lapply(per_sample_loops,
                             function(d) d[d$chrA == d$chrB, , drop = FALSE])

  # ---- Bin both anchors at `resolution`, build consensus ----
  bin_pairs <- vector("list", length(per_sample_loops))
  for (i in seq_along(per_sample_loops)) {
    d <- per_sample_loops[[i]]
    if (nrow(d) == 0L) {
      bin_pairs[[i]] <- data.frame(chr = character(0),
                                   binA = integer(0), binB = integer(0),
                                   stringsAsFactors = FALSE)
      next
    }
    bA <- .bin_anchor(d$startA, d$endA, resolution)
    bB <- .bin_anchor(d$startB, d$endB, resolution)
    swap <- bA > bB
    a_lo <- ifelse(swap, bB, bA)
    a_hi <- ifelse(swap, bA, bB)
    bin_pairs[[i]] <- unique(data.frame(
      chr  = d$chrA,
      binA = a_lo,
      binB = a_hi,
      stringsAsFactors = FALSE
    ))
  }

  all_pairs <- do.call(rbind, bin_pairs)
  all_pairs$key <- paste(all_pairs$chr, all_pairs$binA, all_pairs$binB, sep = ":")
  pair_counts <- table(all_pairs$key)
  keep_keys <- names(pair_counts)[pair_counts >= min_samples]
  if (length(keep_keys) == 0L)
    stop("No consensus loops survived min_samples = ", min_samples)

  uniq <- unique(all_pairs[all_pairs$key %in% keep_keys,
                            c("chr","binA","binB","key"), drop = FALSE])
  ord <- order(uniq$chr, uniq$binA, uniq$binB)
  uniq <- uniq[ord, , drop = FALSE]
  rownames(uniq) <- NULL

  loop_id <- sprintf("L%07d_%s_%d_%d",
                     seq_len(nrow(uniq)), uniq$chr, uniq$binA, uniq$binB)

  anchorA <- GenomicRanges::GRanges(
    seqnames = uniq$chr,
    ranges   = IRanges::IRanges(start = uniq$binA + 1L,
                                end   = uniq$binA + resolution)
  )
  anchorB <- GenomicRanges::GRanges(
    seqnames = uniq$chr,
    ranges   = IRanges::IRanges(start = uniq$binB + 1L,
                                end   = uniq$binB + resolution)
  )
  gi <- InteractionSet::GInteractions(anchorA, anchorB)
  S4Vectors::mcols(gi)$loop_id <- loop_id

  # ---- cis: user-supplied if provided, else auto-computed from .hic ----
  cis_chroms <- unique(uniq$chr)   # chroms with consensus loops; defines the
                                   # set of chroms cis is summed over (consistent
                                   # with what normalize will process)

  cis_vec <- if (!is.null(user_cis)) {
    if (anyNA(user_cis) || any(user_cis <= 0))
      stop("`samples$cis` must be positive numeric with no NAs")
    setNames(as.numeric(user_cis), sample_ids)
  } else {
    message("[", format(Sys.time()), "] computing cis depth from .hic files ",
            "(no `cis` column in samplesheet; ", length(cis_chroms),
            " chroms x ", length(sample_ids), " samples) ...")
    .compute_cis_from_hic(
      hic_paths   = hic_paths,
      sample_ids  = sample_ids,
      chroms      = cis_chroms,
      resolution  = resolution,
      chrom_norms = chrom_norm_per_sample,
      BPPARAM     = BPPARAM
    )
  }

  cis_source_tag <- if (!is.null(user_cis)) "user_provided" else "computed_from_hic"

  meta <- list(
    resolution         = resolution,
    min_samples        = min_samples,
    n_trans_dropped    = sum(n_trans_per_sample),
    n_trans_per_sample = setNames(n_trans_per_sample, sample_ids),
    cis_source         = cis_source_tag,
    cis_chroms_in      = if (length(cis_vec) > 0L) cis_chroms else character(0),
    normalization      = "NONE",
    hic_paths          = setNames(hic_paths, sample_ids)
  )

  obj <- new("loopQTL",
             samples    = sample_ids,
             loops      = gi,
             contacts   = matrix(numeric(0), 0, 0),
             cis        = cis_vec,
             vcf_path   = vcf,
             resolution = resolution,
             phenotype  = matrix(numeric(0), 0, 0),
             meta       = meta)
  methods::validObject(obj)
  obj
}

# Internal: compute per-sample cis depth by summing all straw counts across
# the given chroms. Returns a named numeric vector (one entry per sample).
# Parallel unit of work: one (sample, chrom) straw pull + sum.
.compute_cis_from_hic <- function(hic_paths, sample_ids, chroms,
                                   resolution, chrom_norms,
                                   BPPARAM = BiocParallel::SerialParam()) {
  grid <- expand.grid(sample_i = seq_along(sample_ids),
                      chrom    = chroms,
                      stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS = FALSE)

  one <- function(k) {
    si  <- grid$sample_i[k]
    chr <- grid$chrom[k]
    chr_in_hic <- chrom_norms[[si]](chr)
    dat <- tryCatch(
      strawr::straw("NONE", hic_paths[si], chr_in_hic, chr_in_hic,
                    "BP", resolution),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0L) return(0)
    sum(as.numeric(dat$counts))
  }

  results <- BiocParallel::bplapply(seq_len(nrow(grid)), one, BPPARAM = BPPARAM)
  cis_vec <- setNames(numeric(length(sample_ids)), sample_ids)
  for (k in seq_len(nrow(grid))) {
    cis_vec[grid$sample_i[k]] <- cis_vec[grid$sample_i[k]] + results[[k]]
  }
  cis_vec
}
