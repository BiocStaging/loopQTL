#' Anchor-restricted loop-by-SNP association tests
#'
#' For each consensus loop, selects SNPs overlapping either anchor (within an
#' optional `window` bp flank), reads dosages from the VCF in a single targeted
#' tabix pass, and runs `lm(phenotype ~ dosage)` per (loop, SNP). No covariates
#' or transforms applied here — those are owned by [loopQTL_normalize()].
#'
#' @param obj A [loopQTL-class] after [loopQTL_normalize()].
#' @param window Bp flank added to each anchor bin when selecting SNPs.
#'   `window = 0` keeps SNPs inside the anchor bin only. Larger windows erode
#'   the anchor-restriction rationale; not enforced.
#'
#' @return data.frame with one row per tested (loop, SNP) pair:
#'   `loop_id`, `snp_id`, `chrom`, `pos`, `ref`, `alt`, `anchor` (A/B/both),
#'   `beta`, `se`, `t`, `p`, `n_used`, `maf`.
#'
#' @export
run_QTL <- function(obj, window = 0L) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    length(obj@phenotype) > 0L,
    is.numeric(window), length(window) == 1L, window >= 0
  )
  window <- as.integer(window)

  loops_gi <- obj@loops
  resolution <- obj@resolution
  phenotype <- obj@phenotype
  samples_vec <- obj@samples
  loop_ids <- S4Vectors::mcols(loops_gi)$loop_id

  anchorA <- InteractionSet::anchors(loops_gi, type = "first")
  anchorB <- InteractionSet::anchors(loops_gi, type = "second")

  # ---- Build the read-target ranges (union of A and B, padded by window) ----
  pad <- function(gr) {
    GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr),
      ranges   = IRanges::IRanges(
        start = pmax(1L, GenomicRanges::start(gr) - window),
        end   = GenomicRanges::end(gr) + window
      )
    )
  }
  rangesA <- pad(anchorA)
  rangesB <- pad(anchorB)
  read_targets <- GenomicRanges::reduce(c(rangesA, rangesB))

  # ---- Detect VCF chromosome naming style and rewrite if needed ----
  tab <- Rsamtools::TabixFile(obj@vcf_path)
  vcf_chroms <- Rsamtools::headerTabix(tab)$seqnames
  if (length(vcf_chroms) == 0L)
    stop("Could not read chromosome names from VCF tabix index")
  vcf_norm <- .make_chrom_normalizer(vcf_chroms)
  rt_chroms <- as.character(GenomicRanges::seqnames(read_targets))
  rt_chroms_vcf <- vcf_norm(rt_chroms)
  GenomicRanges::seqlevels(read_targets) <- unique(rt_chroms_vcf)
  GenomicRanges::seqnames(read_targets) <-
    S4Vectors::Rle(rt_chroms_vcf)

  read_targets <- read_targets[as.character(GenomicRanges::seqnames(read_targets)) %in% vcf_chroms]
  if (length(read_targets) == 0L) {
    return(.empty_qtl_result())
  }

  # ---- Single targeted VCF read ----
  param <- VariantAnnotation::ScanVcfParam(
    which = read_targets,
    geno  = "GT",
    info  = NA_character_
  )
  vcf <- VariantAnnotation::readVcf(obj@vcf_path, genome = "unknown", param = param)
  if (length(vcf) == 0L) {
    return(.empty_qtl_result())
  }

  # Reorder VCF samples to the canonical order; error if any missing
  vcf_sample_names <- colnames(VariantAnnotation::geno(vcf)$GT)
  m <- match(samples_vec, vcf_sample_names)
  if (anyNA(m))
    stop("VCF missing samples that were validated at setup: ",
         paste(samples_vec[is.na(m)], collapse = ", "))

  gt <- VariantAnnotation::geno(vcf)$GT[, m, drop = FALSE]
  dosage <- .gt_to_dosage(gt)
  colnames(dosage) <- samples_vec

  rr <- SummarizedExperiment::rowRanges(vcf)
  # also normalize snp chrom names back to user-style for the result
  snp_chrom <- as.character(GenomicRanges::seqnames(rr))
  snp_pos   <- GenomicRanges::start(rr)
  snp_ref   <- as.character(VariantAnnotation::ref(vcf))
  alt_list  <- VariantAnnotation::alt(vcf)
  snp_alt   <- vapply(alt_list,
                      function(a) paste(as.character(a), collapse = ","),
                      character(1L))
  snp_id    <- rownames(vcf)
  if (is.null(snp_id) || any(snp_id == "" | is.na(snp_id)))
    snp_id <- sprintf("%s:%d_%s/%s", snp_chrom, snp_pos, snp_ref, snp_alt)

  # ---- Map SNPs to loops via anchor overlap ----
  # Build matched-naming anchor GRanges for overlap with the VCF rowRanges
  anchorA_vcf <- anchorA
  anchorB_vcf <- anchorB
  GenomicRanges::seqlevels(anchorA_vcf) <- unique(vcf_norm(GenomicRanges::seqlevels(anchorA_vcf)))
  GenomicRanges::seqnames(anchorA_vcf)  <-
    S4Vectors::Rle(vcf_norm(as.character(GenomicRanges::seqnames(anchorA))))
  GenomicRanges::seqlevels(anchorB_vcf) <- unique(vcf_norm(GenomicRanges::seqlevels(anchorB_vcf)))
  GenomicRanges::seqnames(anchorB_vcf)  <-
    S4Vectors::Rle(vcf_norm(as.character(GenomicRanges::seqnames(anchorB))))

  padded_A <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(anchorA_vcf),
    ranges   = IRanges::IRanges(
      start = pmax(1L, GenomicRanges::start(anchorA_vcf) - window),
      end   = GenomicRanges::end(anchorA_vcf) + window
    )
  )
  padded_B <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(anchorB_vcf),
    ranges   = IRanges::IRanges(
      start = pmax(1L, GenomicRanges::start(anchorB_vcf) - window),
      end   = GenomicRanges::end(anchorB_vcf) + window
    )
  )

  hitsA <- GenomicRanges::findOverlaps(padded_A, rr)
  hitsB <- GenomicRanges::findOverlaps(padded_B, rr)
  # combine: per (loop, snp) record which anchor(s) hit
  keyA <- paste(S4Vectors::queryHits(hitsA), S4Vectors::subjectHits(hitsA), sep = ":")
  keyB <- paste(S4Vectors::queryHits(hitsB), S4Vectors::subjectHits(hitsB), sep = ":")
  in_A <- union(keyA, keyB)
  if (length(in_A) == 0L) {
    return(.empty_qtl_result())
  }
  pair <- do.call(rbind, strsplit(in_A, ":", fixed = TRUE))
  loop_i <- as.integer(pair[, 1L])
  snp_j  <- as.integer(pair[, 2L])
  in_setA <- in_A %in% keyA
  in_setB <- in_A %in% keyB
  anchor_lab <- ifelse(in_setA & in_setB, "both",
                       ifelse(in_setA, "A", "B"))

  # ---- Regress each pair ----
  n_pair <- length(loop_i)
  beta <- se <- tstat <- pval <- maf <- rep(NA_real_, n_pair)
  n_used <- rep(NA_integer_, n_pair)
  n_samples <- length(samples_vec)

  for (k in seq_len(n_pair)) {
    li <- loop_i[k]
    sj <- snp_j[k]
    y <- phenotype[li, ]
    x <- dosage[sj, ]
    ok <- is.finite(y) & is.finite(x)
    n_used[k] <- sum(ok)
    if (n_used[k] < 3L) next
    yk <- y[ok]; xk <- x[ok]
    af <- mean(xk) / 2
    maf[k] <- min(af, 1 - af)
    if (stats::sd(xk) == 0) next
    if (stats::sd(yk) == 0) next
    fit <- stats::lm.fit(cbind(1, xk), yk)
    rdf <- fit$df.residual
    if (rdf < 1L) next
    rss <- sum(fit$residuals^2)
    sig2 <- rss / rdf
    # se for slope: sqrt(sig2 * (X'X)^{-1}[2,2])
    xc <- xk - mean(xk)
    sxx <- sum(xc^2)
    if (sxx == 0) next
    b <- fit$coefficients[2L]
    s <- sqrt(sig2 / sxx)
    t_  <- b / s
    p_  <- 2 * stats::pt(-abs(t_), df = rdf)
    beta[k] <- b; se[k] <- s; tstat[k] <- t_; pval[k] <- p_
  }

  data.frame(
    loop_id = loop_ids[loop_i],
    snp_id  = snp_id[snp_j],
    chrom   = snp_chrom[snp_j],
    pos     = snp_pos[snp_j],
    ref     = snp_ref[snp_j],
    alt     = snp_alt[snp_j],
    anchor  = anchor_lab,
    beta    = beta,
    se      = se,
    t       = tstat,
    p       = pval,
    n_used  = n_used,
    maf     = maf,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ---- internal helpers used only by run_QTL ----

.empty_qtl_result <- function() {
  data.frame(
    loop_id = character(0), snp_id = character(0),
    chrom = character(0), pos = integer(0),
    ref = character(0), alt = character(0),
    anchor = character(0),
    beta = numeric(0), se = numeric(0), t = numeric(0), p = numeric(0),
    n_used = integer(0), maf = numeric(0),
    stringsAsFactors = FALSE
  )
}

# Convert a (snps x samples) GT character matrix to numeric dosage (0/1/2).
# Anything not unambiguously diploid biallelic becomes NA.
.gt_to_dosage <- function(gt) {
  # Strip phase, then split on / so "0|1" and "0/1" parse identically.
  raw <- gsub("|", "/", gt, fixed = TRUE)
  raw[raw %in% c("", ".", "./.")] <- NA_character_
  out <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt))
  ok <- !is.na(raw)
  if (!any(ok)) return(out)
  sp <- strsplit(raw[ok], "/", fixed = TRUE)
  bad <- vapply(sp, function(x) length(x) != 2L, logical(1L))
  a1 <- suppressWarnings(as.integer(vapply(sp, `[`, character(1L), 1L)))
  a2 <- suppressWarnings(as.integer(vapply(sp, `[`, character(1L), 2L)))
  bad <- bad | is.na(a1) | is.na(a2) | a1 > 1L | a2 > 1L
  vals <- ifelse(bad, NA_real_, a1 + a2)
  out[ok] <- vals
  out
}
