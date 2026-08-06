#' Extract per-sample data for one or more SNP-loop pairs
#'
#' Returns a long-format data.frame with one row per `(sample, loop_id,
#' snp_id)` combination, carrying the per-sample genotype dosage, the
#' phenotype value at that loop, and the per-sample cis depth.
#'
#' Requires `obj@results` to be populated (i.e., `runQTL()` has run), since
#' SNP positions are looked up from the result table to fetch genotypes.
#'
#' @param obj A [loopQTL-class] with `runQTL()` results in `obj@results`.
#' @param loop_id Character vector of loop_ids to extract.
#' @param snp_id Character vector of snp_ids to extract.
#' @param normal Apply RINT per loop before returning phenotype values
#'   (default `TRUE`, matching `runQTL()`'s default - the values plotted
#'   then match what the regression actually used).
#'
#' @return data.frame with columns `sample`, `loop_id`, `snp_id`, `genotype`
#'   (0/1/2 dosage; NA for missing), `phenotype`, `cis`.
#'
#' @export
loopQTLpairData <- function(obj, loop_id, snp_id, normal = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(loop_id), length(loop_id) > 0L,
        is.character(snp_id), length(snp_id) > 0L,
        is.logical(normal), length(normal) == 1L
    )
    if (nrow(obj@results) == 0L) {
        stop("obj@results is empty. Run runQTL(obj) first to populate it.")
    }

    loop_ids_all <- S4Vectors::mcols(obj@loops)$loop_id
    miss_loops <- setdiff(loop_id, loop_ids_all)
    if (length(miss_loops) > 0L) {
        stop("loop_id(s) not in obj@loops: ", paste(miss_loops, collapse = ", "))
    }

    miss_snps <- setdiff(snp_id, obj@results$snp_id)
    if (length(miss_snps) > 0L) {
        stop("snp_id(s) not in obj@results: ", paste(miss_snps, collapse = ", "))
    }

    # SNP info from obj@results (positions + ref/alt)
    snp_info <- unique(obj@results[
        obj@results$snp_id %in% snp_id,
        c("snp_id", "chrom", "pos", "ref", "alt"),
        drop = FALSE
    ])
    snp_info <- snp_info[match(snp_id, snp_info$snp_id), , drop = FALSE]

    # Genotypes via targeted tabix read  # Normalize chrom naming to whatever the VCF uses.
    tab <- Rsamtools::TabixFile(obj@vcf_path)
    vcf_chroms <- Rsamtools::headerTabix(tab)$seqnames
    vcf_norm <- .make_chrom_normalizer(vcf_chroms)
    snp_chrom_vcf <- vcf_norm(snp_info$chrom)
    which_gr <- GenomicRanges::GRanges(
        seqnames = snp_chrom_vcf,
        ranges   = IRanges::IRanges(start = snp_info$pos, end = snp_info$pos)
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = which_gr, geno = "GT", info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(obj@vcf_path, genome = "unknown", param = param)

    rr <- SummarizedExperiment::rowRanges(vcf)
    vcf_chrom <- as.character(GenomicRanges::seqnames(rr))
    vcf_pos <- GenomicRanges::start(rr)
    vcf_ref <- as.character(VariantAnnotation::ref(vcf))

    # Match by (chrom, pos, ref) to handle multi-allelic loci
    vcf_key <- paste(vcf_chrom, vcf_pos, vcf_ref, sep = ":")
    snp_key <- paste(snp_chrom_vcf, snp_info$pos, snp_info$ref, sep = ":")
    m_vcf <- match(snp_key, vcf_key)
    if (anyNA(m_vcf)) {
        stop(
            "Could not find SNP(s) in VCF at expected (chrom, pos, ref): ",
            paste(snp_id[is.na(m_vcf)], collapse = ", ")
        )
    }

    gt <- VariantAnnotation::geno(vcf)$GT[m_vcf, , drop = FALSE]

    # Reorder VCF samples to canonical
    sm <- match(obj@samples, colnames(gt))
    if (anyNA(sm)) {
        stop(
            "VCF missing samples: ",
            paste(obj@samples[is.na(sm)], collapse = ", ")
        )
    }
    gt <- gt[, sm, drop = FALSE]
    dosage <- .gt_to_dosage(gt)
    rownames(dosage) <- snp_id
    colnames(dosage) <- obj@samples

    # Phenotype matrix for the requested loops
    loop_idx <- match(loop_id, loop_ids_all)
    pheno_full <- if (normal) .rint_matrix(obj@phenotype) else obj@phenotype
    pheno_sub <- pheno_full[loop_idx, , drop = FALSE]
    rownames(pheno_sub) <- loop_id

    # Build long table: n_samples x n_loops x n_snps rows
    n_samples <- length(obj@samples)
    grid <- expand.grid(
        sample_idx = seq_len(n_samples),
        loop_idx = seq_along(loop_id),
        snp_idx = seq_along(snp_id),
        KEEP.OUT.ATTRS = FALSE
    )
    cis_vec <- if (length(obj@cis) > 0L) obj@cis[obj@samples] else rep(NA_real_, n_samples)

    out <- data.frame(
        sample = obj@samples[grid$sample_idx],
        loop_id = loop_id[grid$loop_idx],
        snp_id = snp_id[grid$snp_idx],
        genotype = dosage[cbind(grid$snp_idx, grid$sample_idx)],
        phenotype = pheno_sub[cbind(grid$loop_idx, grid$sample_idx)],
        cis = cis_vec[grid$sample_idx],
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}
