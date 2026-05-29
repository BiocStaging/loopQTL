#' Set up a loopQTL analysis
#'
#' Builds the bin-matched consensus loop set, validates per-sample `.hic`
#' metadata and the VCF header, and records the per-sample `.hic` paths
#' needed by [loopQTL_normalize()]. No contact extraction is performed here
#' (that's the whole-genome pull inside `loopQTL_normalize()`, deliberately
#' kept out of setup so this step is fast and you can iterate on the
#' samplesheet without paying for straw extraction).
#'
#' @param samples data.frame with columns `Sample`, `hic`, `loops`. Optional
#'   `cis` column: if present, those values become `obj@cis` (used downstream
#'   only as a stored QC metric — not by `loopQTL_normalize()`). Absent: `cis`
#'   slot stays empty.
#' @param vcf Path to a bgzipped + tabix-indexed multi-sample VCF.
#' @param resolution Analysis bin size in bp (must exist in every .hic).
#' @param min_samples Keep consensus loops appearing in at least this many
#'   samples (default 2).
#'
#' @return A [loopQTL-class] object with `samples`, `loops`, `vcf_path`,
#'   `resolution`, optionally `cis`, and `meta` populated. `contacts` and
#'   `phenotype` are empty placeholders until normalize fills `phenotype`.
#'
#' @export
loopQTL_setup <- function(samples,
                          vcf,
                          resolution,
                          min_samples = 2L) {

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
  for (i in seq_along(hic_paths)) {
    bp_res <- strawr::readHicBpResolutions(hic_paths[i])
    norms  <- strawr::readHicNormTypes(hic_paths[i])
    if (!(resolution %in% bp_res))
      stop("Resolution ", resolution, " not available in ", hic_paths[i],
           " (available: ", paste(bp_res, collapse = ", "), ")")
    if (!("NONE" %in% norms))
      stop("Normalization 'NONE' not available in ", hic_paths[i])
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

  # ---- cis (user-supplied only; no straw computation in setup) ----
  cis_vec <- if (!is.null(user_cis)) {
    if (anyNA(user_cis) || any(user_cis <= 0))
      stop("`samples$cis` must be positive numeric with no NAs")
    setNames(as.numeric(user_cis), sample_ids)
  } else {
    numeric(0)
  }

  meta <- list(
    resolution         = resolution,
    min_samples        = min_samples,
    n_trans_dropped    = sum(n_trans_per_sample),
    n_trans_per_sample = setNames(n_trans_per_sample, sample_ids),
    cis_source         = if (is.null(user_cis)) "not_provided" else "user_provided",
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
