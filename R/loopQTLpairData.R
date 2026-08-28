#' Extract per-sample data for one or more SNP-loop pairs
#'
#' Returns a long-format data.frame with one row per `(sample, loopID,
#' snpId)` combination, carrying the per-sample genotype dosage, the
#' phenotype value at that loop, and the per-sample cis depth.
#'
#' Requires `obj@results` to be populated (i.e., `runQTL()` has run), since
#' SNP positions are looked up from the result table to fetch genotypes.
#'
#' @param obj A [loopQTL-class] with `runQTL()` results in `obj@results`.
#' @param loopID Character vector of loopIDs to extract.
#' @param snpId Character vector of snpIds to extract.
#' @param normal Apply RINT per loop before returning phenotype values
#'   (default `TRUE`, matching `runQTL()`'s default - the values plotted
#'   then match what the regression actually used).
#'
#' @return data.frame with columns `sample`, `loopID`, `snpId`, `genotype`
#'   (0/1/2 dosage; NA for missing), `phenotype`, `cis`.
#'
#' @export
loopQTLpairData <- function(obj, loopID, snpId, normal = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(loopID), length(loopID) > 0,
        is.character(snpId), length(snpId) > 0,
        is.logical(normal), length(normal) == 1
    )
    if (nrow(obj@results) == 0) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first to populate it.")
    }

    loop_ids_all <- S4Vectors::mcols(obj@loops)$loopID
    missingLoops <- setdiff(loopID, loop_ids_all)
    if (length(missingLoops) > 0) {
        .stopNoCall("loopID(s) not in obj@loops: ", paste(missingLoops, collapse = ", "))
    }

    missingSNPs <- setdiff(snpId, obj@results$snpId)
    if (length(missingSNPs) > 0) {
        .stopNoCall("snpId(s) not in obj@results: ", paste(missingSNPs, collapse = ", "))
    }

    ## SNP info from obj@results (positions + ref/alt)
    snpInfo <- unique(obj@results[
        obj@results$snpId %in% snpId,
        c("snpId", "chrom", "pos", "ref", "alt"),
        drop = FALSE
    ])
    snpInfo <- snpInfo[match(snpId, snpInfo$snpId), , drop = FALSE]

    ## Genotypes via targeted tabix read  # Normalize chrom naming to whatever the VCF uses.
    tab <- Rsamtools::TabixFile(obj@vcfPath)
    vcfChroms <- Rsamtools::headerTabix(tab)$seqnames
    vcfChromNormalizerFunc <- .make_chrom_normalizer(vcfChroms)
    snp_chrom_vcf <- vcfChromNormalizerFunc(snpInfo$chrom)
    whichGr <- GenomicRanges::GRanges(
        seqnames = snp_chrom_vcf,
        ranges   = IRanges::IRanges(start = snpInfo$pos, end = snpInfo$pos)
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = whichGr, geno = "GT", info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(obj@vcfPath, genome = "unknown", param = param)

    rr <- SummarizedExperiment::rowRanges(vcf)
    vcfChrom <- as.character(GenomicRanges::seqnames(rr))
    vcfPos <- GenomicRanges::start(rr)
    vcfRef <- as.character(VariantAnnotation::ref(vcf))

    ## Match by (chrom, pos, ref) to handle multi-allelic loci
    vcfKey <- paste(vcfChrom, vcfPos, vcfRef, sep = ":")
    snpKey <- paste(snp_chrom_vcf, snpInfo$pos, snpInfo$ref, sep = ":")
    snpsInVCF <- match(snpKey, vcfKey)
    if (anyNA(snpsInVCF)) {
        .stopNoCall(
            "Could not find SNP(s) in VCF at expected (chrom, pos, ref): ",
            paste(snpId[is.na(snpsInVCF)], collapse = ", ")
        )
    }

    gt <- VariantAnnotation::geno(vcf)$GT[snpsInVCF, , drop = FALSE]

    ## Reorder VCF samples to canonical
    sm <- match(obj@samples, colnames(gt))
    if (anyNA(sm)) {
        .stopNoCall(
            "VCF missing samples: ",
            paste(obj@samples[is.na(sm)], collapse = ", ")
        )
    }
    gt <- gt[, sm, drop = FALSE]
    dosage <- .gt_to_dosage(gt)
    rownames(dosage) <- snpId
    colnames(dosage) <- obj@samples

    ## Phenotype matrix for the requested loops
    loopIdx <- match(loopID, loop_ids_all)
    phenoFull <- if (normal) .rintMatrix(obj@phenotype) else obj@phenotype
    phenoSub <- phenoFull[loopIdx, , drop = FALSE]
    rownames(phenoSub) <- loopID

    ## Build long table: nSamples x nLoops x n_snps rows
    nSamples <- length(obj@samples)
    grid <- expand.grid(
        sampleIdx = seq_len(nSamples),
        loopIdx = seq_along(loopID),
        snpIdx = seq_along(snpId),
        KEEP.OUT.ATTRS = FALSE
    )
    cisVec <- if (length(obj@cis) > 0) obj@cis[obj@samples] else rep(NA_real_, nSamples)

    out <- data.frame(
        sample = obj@samples[grid$sampleIdx],
        loopID = loopID[grid$loopIdx],
        snpId = snpId[grid$snpIdx],
        genotype = dosage[cbind(grid$snpIdx, grid$sampleIdx)],
        phenotype = phenoSub[cbind(grid$loopIdx, grid$sampleIdx)],
        cis = cisVec[grid$sampleIdx],
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}
