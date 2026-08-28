#' Set up a loopQTL analysis
#'
#' Builds the consensus loop set by **exact-coordinate matching** of
#' per-sample BEDPE calls: a loop is "the same" across samples only if both
#' anchors' (start, end) coordinates match exactly. Counts the number of
#' distinct samples each unique loop appears in, and retains those with
#' support >= `minSamples`. No jitter tolerance, no merging - the
#' previously-used mariner clustering step is gone; experience showed it
#' was washing out significant QTL hits by collapsing distinct loops.
#'
#' Also validates per-sample `.hic` metadata and the VCF header, and
#' records the per-sample `.hic` paths needed by [loopQTLnormalize()].
#'
#' When the samplesheet does not provide a `cis` column, per-sample cis
#' depth is unconditionally auto-computed from the `.hic` files via a
#' one-pass straw aggregation (sums of counts per `(sample, chrom)`, no
#' consensus indexing) - needed by the BNBC logCPM step downstream. If you
#' precompute `cis` per sample yourself (e.g. from a Hi-C QC report) and
#' put it in the samplesheet, that value is used as-is and the .hic pass
#' is skipped.
#'
#' @param samples data.frame with columns `Sample`, `hic`, `loops`. Optional
#'   `cis` column: if present, those values become `obj@cis` directly. If
#'   absent, cis is summed from the `.hic` files here at setup time.
#' @param vcf Path to a bgzipped + tabix-indexed multi-sample VCF.
#' @param resolution Analysis bin size in bp (must exist in every .hic). Used
#'   downstream for straw extraction (normalize) and cis aggregation; NOT
#'   used for consensus matching.
#' @param minSamples Keep consensus loops supported by at least this many
#'   distinct samples (default 2). Support is the count of distinct
#'   samples whose BEDPE file contained that exact (chrA, startA, endA,
#'   chrB, startB, endB) tuple after canonicalising A vs B order.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] that parallelizes the
#'   cis-computing straw pulls. Default [BiocParallel::SerialParam()].
#'
#' @return A [loopQTL-class] object with `samples`, `loops`, `vcfPath`,
#'   `resolution`, `cis` (if provided or computed), and `meta` populated.
#'   `phenotype` is an empty placeholder until normalize fills it.
#'
#' @export
loopQTLsetup <- function(samples,
                         vcf,
                         resolution,
                         minSamples = 2,
                         BPPARAM = BiocParallel::SerialParam()) {
    stopifnot(
        is.data.frame(samples),
        all(c("Sample", "hic", "loops") %in% names(samples)),
        is.character(vcf), length(vcf) == 1,
        is.numeric(resolution), length(resolution) == 1, resolution > 0,
        is.numeric(minSamples), length(minSamples) == 1, minSamples >= 1
    )
    resolution <- as.integer(resolution)
    minSamples <- as.integer(minSamples)

    sampleIds <- as.character(samples$Sample)
    if (anyDuplicated(sampleIds)) {
        .stopNoCall("Sample IDs must be unique")
    }
    hicPaths <- as.character(samples$hic)
    loopPaths <- as.character(samples$loops)
    userCis <- if ("cis" %in% names(samples)) as.numeric(samples$cis) else NULL

    for (p in c(hicPaths, loopPaths)) {
        if (!file.exists(p)) .stopNoCall("File not found: ", p)
    }

    ## VCF validation (header only)
    if (!file.exists(vcf)) .stopNoCall("VCF not found: ", vcf)
    tbi <- paste0(vcf, ".tbi")
    if (!file.exists(tbi)) {
        .stopNoCall(
            "Tabix index not found: ", tbi,
            "\nPlease bgzip + tabix-index your VCF (e.g. `tabix -p vcf my.vcf.gz`)."
        )
    }

    vcfSamples <- VariantAnnotation::samples(
        VariantAnnotation::scanVcfHeader(vcf)
    )
    missing_in_vcf <- setdiff(sampleIds, vcfSamples)
    if (length(missing_in_vcf) > 0) {
        message(
            "Dropping ", length(missing_in_vcf),
            " samplesheet sample(s) not in VCF: ",
            paste(missing_in_vcf, collapse = ", ")
        )
        keep <- !(sampleIds %in% missing_in_vcf)
        sampleIds <- sampleIds[keep]
        hicPaths <- hicPaths[keep]
        loopPaths <- loopPaths[keep]
        if (!is.null(userCis)) userCis <- userCis[keep]
    }
    if (length(sampleIds) < 2) {
        .stopNoCall(
            "Fewer than 2 samples overlap between samplesheet and VCF; ",
            "cannot build consensus / run QTL."
        )
    }
    message(
        "Proceeding with ", length(sampleIds),
        " sample(s) present in both samplesheet and VCF."
    )

    ## Per-.hic metadata + resolution / NONE check (header reads only)  # Also build per-sample chromosome-naming normalizers so we can map a
    ## user-style chrom string (e.g. "1") to whatever this .hic file uses
    ## (e.g. "chr1") when summing cis.
    chrom_norm_per_sample <- vector("list", length(hicPaths))
    hic_chroms_per_sample <- vector("list", length(hicPaths))
    for (i in seq_along(hicPaths)) {
        chromsI <- strawr::readHicChroms(hicPaths[i])$name
        bpRes <- strawr::readHicBpResolutions(hicPaths[i])
        norms <- strawr::readHicNormTypes(hicPaths[i])
        if (!(resolution %in% bpRes)) {
            .stopNoCall(
                "Resolution ", resolution, " not available in ", hicPaths[i],
                " (available: ", paste(bpRes, collapse = ", "), ")"
            )
        }
        if (!("NONE" %in% norms)) {
            .stopNoCall("Normalization 'NONE' not available in ", hicPaths[i])
        }
        chrom_norm_per_sample[[i]] <- .make_chrom_normalizer(chromsI)
        hic_chroms_per_sample[[i]] <- chromsI
    }

    ## Read and partition loops
    per_sample_loops <- lapply(loopPaths, .read_loops_bedpe)
    n_trans_per_sample <- vapply(
        per_sample_loops,
        function(d) sum(d$chrA != d$chrB),
        integer(1)
    )
    per_sample_loops <- lapply(
        per_sample_loops,
        function(d) d[d$chrA == d$chrB, , drop = FALSE]
    )

    ## Exact-coordinate consensus  # Canonicalise (smaller anchor first), build a string key per loop, then
    ## count distinct samples per key. Keep keys with support >= minSamples.
    ## No clustering, no jitter tolerance: two loops are the same only if
    ## (chrA, startA, endA, chrB, startB, endB) match exactly.
    canonicalise <- function(d) {
        if (nrow(d) == 0) {
            return(d)
        }
        swap <- d$startA > d$startB
        if (any(swap)) {
            ts <- d$startA[swap]
            te <- d$endA[swap]
            d$startA[swap] <- d$startB[swap]
            d$endA[swap] <- d$endB[swap]
            d$startB[swap] <- ts
            d$endB[swap] <- te
        }
        d
    }
    per_sample_loops <- lapply(per_sample_loops, canonicalise)
    names(per_sample_loops) <- sampleIds

    per_sample_n <- vapply(per_sample_loops, nrow, integer(1))
    if (sum(per_sample_n > 0) < 2) {
        .stopNoCall("Need at least 2 samples with at least one intra-chrom loop call")
    }
    if (any(per_sample_n == 0)) {
        message(
            "Dropping ", sum(per_sample_n == 0),
            " sample(s) with no usable loops: ",
            paste(sampleIds[per_sample_n == 0], collapse = ", ")
        )
    }

    message(
        "exact-coord consensus across ",
        sum(per_sample_n > 0), " samples ..."
    )
    ## Build a long data.frame: one row per per-sample call, with key + sample
    callChunks <- lapply(seq_along(per_sample_loops), function(si) {
        d <- per_sample_loops[[si]]
        if (nrow(d) == 0) {
            return(NULL)
        }
        data.frame(
            key = paste(d$chrA, d$startA, d$endA, d$startB, d$endB, sep = ":"),
            chrA = d$chrA,
            startA = d$startA, endA = d$endA,
            startB = d$startB, endB = d$endB,
            sample = sampleIds[si],
            stringsAsFactors = FALSE
        )
    })
    allCalls <- do.call(rbind, callChunks)
    ## Distinct-sample support per key (count unique samples after de-duping
    ## within-sample repeats)
    uniqCall <- unique(allCalls[, c("key", "sample")])
    support <- table(uniqCall$key)
    keepKeys <- names(support)[support >= minSamples]
    message(
        length(support), " unique loop tuples; ",
        length(keepKeys), " survive minSamples = ", minSamples
    )
    if (length(keepKeys) == 0) {
        .stopNoCall("No consensus loops survived minSamples = ", minSamples)
    }

    ## Pick the first occurrence of each kept key for its coords
    firstIdx <- match(keepKeys, allCalls$key)
    cons <- allCalls[firstIdx, c("chrA", "startA", "endA", "startB", "endB"),
        drop = FALSE
    ]
    ## Stable order: chrom, anchorA start, anchorB start
    ord <- order(cons$chrA, cons$startA, cons$startB)
    cons <- cons[ord, , drop = FALSE]
    rownames(cons) <- NULL
    nSupport <- as.integer(support[paste(cons$chrA, cons$startA, cons$endA,
        cons$startB, cons$endB,
        sep = ":"
    )])

    ## Build GInteractions (1-based inclusive ranges)
    gi <- InteractionSet::GInteractions(
        GenomicRanges::GRanges(
            cons$chrA,
            IRanges::IRanges(cons$startA + 1, cons$endA)
        ),
        GenomicRanges::GRanges(
            cons$chrA,
            IRanges::IRanges(cons$startB + 1, cons$endB)
        )
    )

    ## Stable loopID from anchor coordinates (0-based BEDPE-style starts)
    loopID <- sprintf(
        "L%07d_%s_%d_%d",
        seq_along(gi),
        cons$chrA,
        cons$startA,
        cons$startB
    )
    S4Vectors::mcols(gi)$loopID <- loopID
    S4Vectors::mcols(gi)$nSamples <- nSupport

    repA <- InteractionSet::anchors(gi, type = "first")

    ## cis: user-supplied if provided, else auto-computed from .hic
    cisChroms <- unique(as.character(GenomicRanges::seqnames(repA))) # chroms with consensus loops; defines the
    ## set of chroms cis is summed over (consistent
    ## with what normalize will process)

    cisVec <- if (!is.null(userCis)) {
        if (anyNA(userCis) || any(userCis <= 0)) {
            .stopNoCall("`samples$cis` must be positive numeric with no NAs")
        }
        setNames(as.numeric(userCis), sampleIds)
    } else {
        message(
            "computing cis depth from .hic files ",
            "(no `cis` column in samplesheet; ", length(cisChroms),
            " chroms x ", length(sampleIds), " samples) ..."
        )
        .compute_cis_from_hic(
            hicPaths   = hicPaths,
            sampleIds  = sampleIds,
            chroms      = cisChroms,
            resolution  = resolution,
            chromNorms = chrom_norm_per_sample,
            BPPARAM     = BPPARAM
        )
    }

    cis_source_tag <- if (!is.null(userCis)) "userProvided" else "computed_from_hic"

    meta <- list(
        resolution         = resolution,
        minSamples        = minSamples,
        n_trans_dropped    = sum(n_trans_per_sample),
        n_trans_per_sample = setNames(n_trans_per_sample, sampleIds),
        cisSource         = cis_source_tag,
        cis_chroms_in      = if (length(cisVec) > 0) cisChroms else character(0),
        normalization      = "NONE",
        hicPaths          = setNames(hicPaths, sampleIds)
    )

    obj <- new("loopQTL",
        samples    = sampleIds,
        loops      = gi,
        cis        = cisVec,
        vcfPath   = vcf,
        resolution = resolution,
        phenotype  = matrix(numeric(0), 0, 0),
        meta       = meta
    )
    methods::validObject(obj)
    obj
}

## Internal: compute per-sample cis depth by summing all straw counts across
## the given chroms. Returns a named numeric vector (one entry per sample).
## Parallel unit of work: one (sample, chrom) straw pull + sum.
.compute_cis_from_hic <- function(hicPaths, sampleIds, chroms,
                                  resolution, chromNorms,
                                  BPPARAM = BiocParallel::SerialParam()) {
    grid <- expand.grid(
        sampleI = seq_along(sampleIds),
        chrom = chroms,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    one <- function(k) {
        si <- grid$sampleI[k]
        chr <- grid$chrom[k]
        chr_in_hic <- chromNorms[[si]](chr)
        dat <- tryCatch(
            strawr::straw(
                "NONE", hicPaths[si], chr_in_hic, chr_in_hic,
                "BP", resolution
            ),
            error = function(e) NULL
        )
        if (is.null(dat) || nrow(dat) == 0) {
            return(0)
        }
        sum(as.numeric(dat$counts))
    }

    results <- BiocParallel::bplapply(seq_len(nrow(grid)), one, BPPARAM = BPPARAM)
    cisVec <- setNames(numeric(length(sampleIds)), sampleIds)
    for (k in seq_len(nrow(grid))) {
        cisVec[grid$sampleI[k]] <- cisVec[grid$sampleI[k]] + results[[k]]
    }
    cisVec
}
