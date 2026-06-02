#' Anchor-restricted loop-by-SNP association tests via MatrixEQTL
#'
#' Reads anchor-window SNPs from the VCF in a single targeted tabix pass, then
#' runs MatrixEQTL twice — once with anchor-A as the gene position and once
#' with anchor-B. Results are outer-joined on `(loop_id, snp_id)` and labeled
#' `anchor = A / B / both`.
#'
#' RINT (Blom-offset rank-inverse-normal transform) is applied per loop just
#' before the regression when `normal = TRUE` (default). Zero-variance loops
#' become all-NA phenotype rows and are skipped by MatrixEQTL. The raw,
#' pre-RINT BNBC/fastlo output in `obj@phenotype` is left untouched, so you
#' can A/B-compare the two by toggling this flag.
#'
#' Covariates: if `covariates = NULL`, K phenotype PCs are auto-chosen via
#' PCAForQTL (default `K_method = "elbow"`). PCA is run on the same phenotype
#' form (RINT'd or not) the regression will use. Otherwise the user-supplied
#' matrix is used as-is.
#'
#' @param obj A [loopQTL-class] after [loopQTL_normalize()].
#' @param window Bp flank added to each anchor bin when matching SNPs. Default 0.
#' @param covariates Optional `samples x k` matrix/data.frame. If `NULL`,
#'   auto-chosen via PCAForQTL.
#' @param normal Apply RINT per loop before regression. Default `TRUE`.
#' @param K_method `"elbow"` or `"BE"`; passed to [loopQTL_choose_K()] when
#'   covariates is `NULL`. Default `"elbow"`.
#' @param pvOutputThreshold MatrixEQTL p-value threshold. Default 1 (return
#'   every pair). Lower it to shrink the output for very large runs.
#'
#' @return data.frame with one row per tested (loop, SNP) pair: `loop_id`,
#'   `snp_id`, `chrom`, `pos`, `ref`, `alt`, `anchor`, `beta`, `se`, `t`, `p`,
#'   `FDR_cis`, `n_used`, `maf`.
#'
#' @export
run_QTL <- function(obj,
                    window            = 0L,
                    covariates        = NULL,
                    normal            = TRUE,
                    K_method          = c("elbow", "BE"),
                    pvOutputThreshold = 1) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    length(obj@phenotype) > 0L,
    is.numeric(window), length(window) == 1L, window >= 0,
    is.logical(normal), length(normal) == 1L, !is.na(normal)
  )
  window <- as.integer(window)
  K_method <- match.arg(K_method)

  loops_gi    <- obj@loops
  samples_vec <- obj@samples
  loop_ids    <- S4Vectors::mcols(loops_gi)$loop_id
  anchorA <- InteractionSet::anchors(loops_gi, type = "first")
  anchorB <- InteractionSet::anchors(loops_gi, type = "second")

  # ---- Apply RINT if requested. obj@phenotype is left untouched (the local
  # `phenotype` variable is what flows into MatrixEQTL and PCA). ----
  phenotype <- if (normal) .rint_matrix(obj@phenotype) else obj@phenotype
  n_zero_var <- sum(rowSums(!is.na(phenotype)) == 0L) -
                sum(rowSums(!is.na(obj@phenotype)) == 0L)
  # n_zero_var = loops newly NA'd by RINT (had data but constant after RINT).
  # Originally-NA loops aren't counted here.

  # ---- Build union of (anchor +- window) target ranges for the VCF read ----
  pad <- function(gr) {
    GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(gr),
      ranges   = IRanges::IRanges(
        start = pmax(1L, GenomicRanges::start(gr) - window),
        end   = GenomicRanges::end(gr) + window
      )
    )
  }
  read_targets <- GenomicRanges::reduce(c(pad(anchorA), pad(anchorB)))

  # Detect VCF chrom naming style and rewrite the read targets to match
  tab <- Rsamtools::TabixFile(obj@vcf_path)
  vcf_chroms <- Rsamtools::headerTabix(tab)$seqnames
  if (length(vcf_chroms) == 0L)
    stop("Could not read chromosome names from VCF tabix index")
  vcf_norm <- .make_chrom_normalizer(vcf_chroms)
  rt_chroms_vcf <- vcf_norm(as.character(GenomicRanges::seqnames(read_targets)))
  keep <- rt_chroms_vcf %in% vcf_chroms
  if (!any(keep)) return(.empty_qtl_result())
  read_targets <- GenomicRanges::GRanges(
    seqnames = rt_chroms_vcf[keep],
    ranges   = IRanges::IRanges(
      start = GenomicRanges::start(read_targets)[keep],
      end   = GenomicRanges::end(read_targets)[keep]
    )
  )

  # ---- Single targeted VCF read ----
  param <- VariantAnnotation::ScanVcfParam(
    which = read_targets, geno = "GT", info = NA_character_
  )
  vcf <- VariantAnnotation::readVcf(obj@vcf_path, genome = "unknown", param = param)
  if (length(vcf) == 0L) return(.empty_qtl_result())

  vcf_sample_names <- colnames(VariantAnnotation::geno(vcf)$GT)
  m <- match(samples_vec, vcf_sample_names)
  if (anyNA(m))
    stop("VCF missing samples that were validated at setup: ",
         paste(samples_vec[is.na(m)], collapse = ", "))

  gt <- VariantAnnotation::geno(vcf)$GT[, m, drop = FALSE]
  dosage <- .gt_to_dosage(gt)
  colnames(dosage) <- samples_vec

  rr <- SummarizedExperiment::rowRanges(vcf)
  snp_chrom_vcf <- as.character(GenomicRanges::seqnames(rr))   # in VCF naming
  snp_pos       <- GenomicRanges::start(rr)
  snp_ref       <- as.character(VariantAnnotation::ref(vcf))
  alt_list      <- VariantAnnotation::alt(vcf)
  snp_alt       <- vapply(alt_list,
                          function(a) paste(as.character(a), collapse = ","),
                          character(1L))
  snp_id_raw <- rownames(vcf)
  if (is.null(snp_id_raw) || any(snp_id_raw == "" | is.na(snp_id_raw)))
    snp_id_raw <- sprintf("%s:%d_%s/%s",
                          snp_chrom_vcf, snp_pos, snp_ref, snp_alt)

  # MAF + drop SNPs with no informative variation
  af <- rowMeans(dosage, na.rm = TRUE) / 2
  maf_all <- pmin(af, 1 - af)
  snp_keep <- is.finite(maf_all) & maf_all > 0
  if (!any(snp_keep)) return(.empty_qtl_result())
  dosage     <- dosage[snp_keep, , drop = FALSE]
  snp_id     <- snp_id_raw[snp_keep]
  snp_chrom_vcf <- snp_chrom_vcf[snp_keep]
  snp_pos       <- snp_pos[snp_keep]
  snp_ref       <- snp_ref[snp_keep]
  snp_alt       <- snp_alt[snp_keep]
  maf_kept      <- maf_all[snp_keep]
  rownames(dosage) <- snp_id

  # ---- Resolve covariates ----
  if (is.null(covariates)) {
    # Use the already-RINT'd-or-not `phenotype` matrix to avoid double-RINT.
    K_res <- .pca_for_qtl(phenotype, samples = samples_vec, K_method = K_method)
    cvrt_mat <- t(K_res$PCs)     # MatrixEQTL wants k x n_samples
    obj_K_meta <- list(K = K_res$K, K_method = K_res$method, K_source = "auto_PCAForQTL")
  } else {
    cov_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
    if (!is.null(rownames(cov_df))) {
      mm <- match(samples_vec, rownames(cov_df))
      if (anyNA(mm))
        stop("`covariates` missing rows for samples: ",
             paste(samples_vec[is.na(mm)], collapse = ", "))
      cov_df <- cov_df[mm, , drop = FALSE]
    } else if (nrow(cov_df) != length(samples_vec)) {
      stop("`covariates` has ", nrow(cov_df), " rows but ",
           length(samples_vec), " samples; provide rownames or matched order")
    }
    cvrt_mat <- t(as.matrix(cov_df))
    colnames(cvrt_mat) <- samples_vec
    obj_K_meta <- list(K = nrow(cvrt_mat), K_method = NA_character_, K_source = "user")
  }
  obj@meta$run_QTL_K        <- obj_K_meta
  obj@meta$run_QTL_normal   <- normal
  obj@meta$run_QTL_n_zero_var <- max(0L, n_zero_var)

  # ---- Build SlicedData inputs for MatrixEQTL ----
  snps_sd <- .as_sliced(dosage)
  gene_sd <- .as_sliced(phenotype, rownames_default = loop_ids)
  cvrt_sd <- .as_sliced(cvrt_mat)

  # SNP positions (use original loop chrom style for output, but MatrixEQTL
  # uses position tables only for cis-distance arithmetic, so we feed it the
  # vcf-style chrom names and matching anchor chrom names)
  loop_chrom_vcf <- vcf_norm(as.character(GenomicRanges::seqnames(anchorA)))
  snpspos <- data.frame(
    snpid = snp_id,
    chr   = snp_chrom_vcf,
    pos   = snp_pos,
    stringsAsFactors = FALSE
  )
  genepos_A <- data.frame(
    geneid = loop_ids,
    chr    = loop_chrom_vcf,
    left   = GenomicRanges::start(anchorA),
    right  = GenomicRanges::end(anchorA),
    stringsAsFactors = FALSE
  )
  genepos_B <- data.frame(
    geneid = loop_ids,
    chr    = loop_chrom_vcf,
    left   = GenomicRanges::start(anchorB),
    right  = GenomicRanges::end(anchorB),
    stringsAsFactors = FALSE
  )

  # ---- Run MatrixEQTL twice ----
  out_A <- .run_matrix_eqtl(snps_sd, gene_sd, cvrt_sd,
                            snpspos, genepos_A,
                            cisDist = window,
                            pvOutputThreshold = pvOutputThreshold)
  out_B <- .run_matrix_eqtl(snps_sd, gene_sd, cvrt_sd,
                            snpspos, genepos_B,
                            cisDist = window,
                            pvOutputThreshold = pvOutputThreshold)

  # Merge on (loop_id = gene, snp_id = snps); label anchor A / B / both.
  merged <- .merge_anchor_runs(out_A, out_B)
  if (nrow(merged) == 0L) return(.empty_qtl_result())

  # Attach SNP metadata + MAF
  snp_meta <- data.frame(
    snp_id = snp_id, chrom = snp_chrom_vcf, pos = snp_pos,
    ref = snp_ref, alt = snp_alt, maf = maf_kept,
    stringsAsFactors = FALSE
  )
  out <- merge(merged, snp_meta, by = "snp_id", all.x = TRUE, sort = FALSE)
  out$n_used <- rowSums(!is.na(dosage[match(out$snp_id, rownames(dosage)), , drop = FALSE]))

  # Final column order
  out <- out[, c("loop_id", "snp_id", "chrom", "pos", "ref", "alt", "anchor",
                 "beta", "se", "t", "p", "FDR_cis", "n_used", "maf")]
  rownames(out) <- NULL
  out
}

# -------------------------------------------------------------------------- #
# Internal helpers
# -------------------------------------------------------------------------- #

.empty_qtl_result <- function() {
  data.frame(
    loop_id = character(0), snp_id = character(0),
    chrom = character(0), pos = integer(0),
    ref = character(0), alt = character(0),
    anchor = character(0),
    beta = numeric(0), se = numeric(0), t = numeric(0), p = numeric(0),
    FDR_cis = numeric(0), n_used = integer(0), maf = numeric(0),
    stringsAsFactors = FALSE
  )
}

# Convert a (snps x samples) GT character matrix to numeric dosage (0/1/2).
# Anything beyond biallelic 0/1 becomes NA.
.gt_to_dosage <- function(gt) {
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

# Wrap a matrix in a MatrixEQTL::SlicedData.
.as_sliced <- function(mat, rownames_default = NULL) {
  if (!is.null(rownames_default) && is.null(rownames(mat)))
    rownames(mat) <- rownames_default
  sd <- MatrixEQTL::SlicedData$new()
  sd$CreateFromMatrix(as.matrix(mat))
  sd
}

# Run MatrixEQTL once and return a tidy data.frame of cis hits.
.run_matrix_eqtl <- function(snps_sd, gene_sd, cvrt_sd,
                             snpspos, genepos,
                             cisDist, pvOutputThreshold) {
  tmp_cis <- tempfile()
  on.exit(unlink(tmp_cis), add = TRUE)
  res <- MatrixEQTL::Matrix_eQTL_main(
    snps                = snps_sd,
    gene                = gene_sd,
    cvrt                = cvrt_sd,
    output_file_name    = NULL,
    pvOutputThreshold   = 0,                 # disable trans output
    useModel            = MatrixEQTL::modelLINEAR,
    errorCovariance     = numeric(0),
    verbose             = FALSE,
    output_file_name.cis = tmp_cis,
    pvOutputThreshold.cis = pvOutputThreshold,
    snpspos             = snpspos,
    genepos             = genepos,
    cisDist             = cisDist,
    pvalue.hist         = FALSE,
    min.pv.by.genesnp   = FALSE,
    noFDRsaveMemory     = FALSE
  )
  cis <- res$cis$eqtls
  if (is.null(cis) || nrow(cis) == 0L)
    return(data.frame(
      loop_id = character(0), snp_id = character(0),
      beta = numeric(0), se = numeric(0), t = numeric(0),
      p = numeric(0), FDR_cis = numeric(0),
      stringsAsFactors = FALSE
    ))
  # MatrixEQTL columns: snps, gene, statistic, pvalue, FDR, beta
  # Some versions also include `se` directly; if not, derive: se = beta / statistic.
  se_col <- if ("se" %in% names(cis)) cis$se else cis$beta / cis$statistic
  data.frame(
    loop_id = as.character(cis$gene),
    snp_id  = as.character(cis$snps),
    beta    = as.numeric(cis$beta),
    se      = as.numeric(se_col),
    t       = as.numeric(cis$statistic),
    p       = as.numeric(cis$pvalue),
    FDR_cis = as.numeric(cis$FDR),
    stringsAsFactors = FALSE
  )
}

# Outer-join two anchor runs on (loop_id, snp_id) and label `anchor`.
# When a pair appears in both runs the model is identical, so we just keep one
# row and label it `both`.
.merge_anchor_runs <- function(out_A, out_B) {
  if (nrow(out_A) == 0L && nrow(out_B) == 0L)
    return(data.frame(
      loop_id = character(0), snp_id = character(0),
      anchor = character(0),
      beta = numeric(0), se = numeric(0), t = numeric(0),
      p = numeric(0), FDR_cis = numeric(0),
      stringsAsFactors = FALSE
    ))
  key <- function(d) paste(d$loop_id, d$snp_id, sep = "::")
  kA <- if (nrow(out_A)) key(out_A) else character(0)
  kB <- if (nrow(out_B)) key(out_B) else character(0)

  in_both <- intersect(kA, kB)
  only_A  <- setdiff(kA, kB)
  only_B  <- setdiff(kB, kA)

  rows <- list()
  if (length(in_both) > 0L) {
    pick <- match(in_both, kA)
    r <- out_A[pick, , drop = FALSE]
    r$anchor <- "both"
    rows[[length(rows) + 1L]] <- r
  }
  if (length(only_A) > 0L) {
    pick <- match(only_A, kA)
    r <- out_A[pick, , drop = FALSE]
    r$anchor <- "A"
    rows[[length(rows) + 1L]] <- r
  }
  if (length(only_B) > 0L) {
    pick <- match(only_B, kB)
    r <- out_B[pick, , drop = FALSE]
    r$anchor <- "B"
    rows[[length(rows) + 1L]] <- r
  }
  do.call(rbind, rows)
}
