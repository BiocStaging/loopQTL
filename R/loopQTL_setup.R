#' Set up a loopQTL analysis
#'
#' Reads per-sample `.hic` files, builds a bin-matched consensus loop set, and
#' assembles a (loops x samples) contact matrix at fixed resolution. Inter-
#' chromosomal loops are dropped (count recorded in `meta`). The VCF header is
#' read here for validation; genotypes are read lazily in [run_QTL()].
#'
#' @param samples data.frame with columns `Sample`, `hic`, `loops`, and
#'   optionally `cis`. If `cis` is absent it is computed from the same straw
#'   pulls used to build the contact matrix.
#' @param vcf Path to a bgzipped + tabix-indexed multi-sample VCF.
#' @param resolution Analysis bin size in bp (must exist in every .hic).
#' @param min_samples Keep consensus loops appearing in at least this many
#'   samples (default 2).
#' @param BPPARAM A [BiocParallel::BiocParallelParam] for parallel
#'   (sample, chromosome) extraction. Default [BiocParallel::SerialParam()].
#'
#' @return A [loopQTL-class] object with `samples`, `loops`, `contacts`, `cis`,
#'   `vcf_path`, `resolution`, and `meta` slots filled. `phenotype` is empty
#'   until [loopQTL_normalize()] runs.
#'
#' @export
loopQTL_setup <- function(samples,
                          vcf,
                          resolution,
                          min_samples = 2L,
                          BPPARAM = BiocParallel::SerialParam()) {

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

  # ---- Per-.hic metadata + resolution / NONE check ----
  hic_chroms_per_sample <- vector("list", length(hic_paths))
  chrom_norm_per_sample <- vector("list", length(hic_paths))
  for (i in seq_along(hic_paths)) {
    chroms <- strawr::readHicChroms(hic_paths[i])$name
    bp_res <- strawr::readHicBpResolutions(hic_paths[i])
    norms  <- strawr::readHicNormTypes(hic_paths[i])
    if (!(resolution %in% bp_res))
      stop("Resolution ", resolution, " not available in ", hic_paths[i],
           " (available: ", paste(bp_res, collapse = ", "), ")")
    if (!("NONE" %in% norms))
      stop("Normalization 'NONE' not available in ", hic_paths[i])
    hic_chroms_per_sample[[i]] <- chroms
    chrom_norm_per_sample[[i]] <- .make_chrom_normalizer(chroms)
  }

  # ---- Read and partition loops ----
  per_sample_loops <- lapply(loop_paths, .read_loops_bedpe)
  n_trans_per_sample <- vapply(per_sample_loops,
                               function(d) sum(d$chrA != d$chrB),
                               integer(1L))
  per_sample_loops <- lapply(per_sample_loops, function(d) d[d$chrA == d$chrB, , drop = FALSE])

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
    # canonicalize so binA <= binB
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
  # stable order: chr, binA, binB
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

  # key -> consensus row index, for fast lookup during extraction
  key_to_idx <- setNames(seq_along(uniq$key), uniq$key)

  # consensus loops by chromosome (for indexing straw output per chrom)
  by_chr <- split(seq_len(nrow(uniq)), uniq$chr)

  # ---- Extract contacts: one straw call per (sample, chromosome) ----
  grid <- expand.grid(sample_i = seq_along(sample_ids),
                      chrom    = names(by_chr),
                      stringsAsFactors = FALSE,
                      KEEP.OUT.ATTRS   = FALSE)

  extract_one <- function(k) {
    si  <- grid$sample_i[k]
    chr <- grid$chrom[k]
    rows <- by_chr[[chr]]
    chr_in_hic <- chrom_norm_per_sample[[si]](chr)
    if (!(chr_in_hic %in% hic_chroms_per_sample[[si]]))
      return(list(rows = rows,
                  counts = rep(0, length(rows)),
                  cis = 0))
    dat <- tryCatch(
      strawr::straw("NONE", hic_paths[si], chr_in_hic, chr_in_hic,
                    "BP", resolution),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0L) {
      return(list(rows = rows,
                  counts = rep(0, length(rows)),
                  cis = 0))
    }
    # canonicalize x <= y for lookup
    swap <- dat$x > dat$y
    if (any(swap)) {
      tmp <- dat$x[swap]; dat$x[swap] <- dat$y[swap]; dat$y[swap] <- tmp
    }
    keys <- paste(chr, dat$x, dat$y, sep = ":")
    want <- uniq$key[rows]
    m <- match(want, keys)
    counts <- ifelse(is.na(m), 0, dat$counts[m])
    list(rows = rows,
         counts = counts,
         cis = sum(dat$counts))
  }

  results <- BiocParallel::bplapply(seq_len(nrow(grid)),
                                    extract_one,
                                    BPPARAM = BPPARAM)

  contacts <- matrix(0, nrow = nrow(uniq), ncol = length(sample_ids),
                     dimnames = list(NULL, sample_ids))
  cis_acc <- setNames(numeric(length(sample_ids)), sample_ids)
  for (k in seq_len(nrow(grid))) {
    si <- grid$sample_i[k]
    contacts[results[[k]]$rows, si] <- results[[k]]$counts
    cis_acc[si] <- cis_acc[si] + results[[k]]$cis
  }

  # ---- cis depth (user-provided overrides computed) ----
  cis_vec <- if (!is.null(user_cis)) {
    if (anyNA(user_cis) || any(user_cis <= 0))
      stop("`samples$cis` must be positive numeric with no NAs")
    setNames(as.numeric(user_cis), sample_ids)
  } else {
    cis_acc
  }

  meta <- list(
    resolution        = resolution,
    min_samples       = min_samples,
    n_trans_dropped   = sum(n_trans_per_sample),
    n_trans_per_sample = setNames(n_trans_per_sample, sample_ids),
    cis_source        = if (is.null(user_cis)) "computed_from_hic" else "user_provided",
    normalization     = "NONE",
    hic_paths         = setNames(hic_paths, sample_ids)   # for loopQTL_normalize's whole-genome straw pull
  )

  obj <- new("loopQTL",
             samples    = sample_ids,
             loops      = gi,
             contacts   = contacts,
             cis        = cis_vec,
             vcf_path   = vcf,
             resolution = resolution,
             phenotype  = matrix(numeric(0), 0, 0),
             meta       = meta)
  methods::validObject(obj)
  obj
}
