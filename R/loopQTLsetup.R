#' Set up a loopQTL analysis
#'
#' Builds the consensus loop set by **exact-coordinate matching** of
#' per-sample BEDPE calls: a loop is "the same" across samples only if both
#' anchors' (start, end) coordinates match exactly. Counts the number of
#' distinct samples each unique loop appears in, and retains those with
#' support >= `min_samples`. No jitter tolerance, no merging - the
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
#' @param min_samples Keep consensus loops supported by at least this many
#'   distinct samples (default 2). Support is the count of distinct
#'   samples whose BEDPE file contained that exact (chrA, startA, endA,
#'   chrB, startB, endB) tuple after canonicalising A vs B order.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] that parallelizes the
#'   cis-computing straw pulls. Default [BiocParallel::SerialParam()].
#'
#' @return A [loopQTL-class] object with `samples`, `loops`, `vcf_path`,
#'   `resolution`, `cis` (if provided or computed), and `meta` populated.
#'   `phenotype` is an empty placeholder until normalize fills it.
#'
#' @export
loopQTLsetup <- function(samples,
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
    resolution <- as.integer(resolution)
    min_samples <- as.integer(min_samples)

    sample_ids <- as.character(samples$Sample)
    if (anyDuplicated(sample_ids)) {
        .stopNoCall("Sample IDs must be unique")
    }
    hic_paths <- as.character(samples$hic)
    loop_paths <- as.character(samples$loops)
    user_cis <- if ("cis" %in% names(samples)) as.numeric(samples$cis) else NULL

    for (p in c(hic_paths, loop_paths)) {
        if (!file.exists(p)) .stopNoCall("File not found: ", p)
    }

    # VCF validation (header only)
    if (!file.exists(vcf)) .stopNoCall("VCF not found: ", vcf)
    tbi <- paste0(vcf, ".tbi")
    if (!file.exists(tbi)) {
        .stopNoCall(
            "Tabix index not found: ", tbi,
            "\nPlease bgzip + tabix-index your VCF (e.g. `tabix -p vcf my.vcf.gz`)."
        )
    }

    vcf_samples <- VariantAnnotation::samples(
        VariantAnnotation::scanVcfHeader(vcf)
    )
    missing_in_vcf <- setdiff(sample_ids, vcf_samples)
    if (length(missing_in_vcf) > 0L) {
        message(
            "Dropping ", length(missing_in_vcf),
            " samplesheet sample(s) not in VCF: ",
            paste(missing_in_vcf, collapse = ", ")
        )
        keep <- !(sample_ids %in% missing_in_vcf)
        sample_ids <- sample_ids[keep]
        hic_paths <- hic_paths[keep]
        loop_paths <- loop_paths[keep]
        if (!is.null(user_cis)) user_cis <- user_cis[keep]
    }
    if (length(sample_ids) < 2L) {
        .stopNoCall(
            "Fewer than 2 samples overlap between samplesheet and VCF; ",
            "cannot build consensus / run QTL."
        )
    }
    message(
        "Proceeding with ", length(sample_ids),
        " sample(s) present in both samplesheet and VCF."
    )

    # Per-.hic metadata + resolution / NONE check (header reads only)  # Also build per-sample chromosome-naming normalizers so we can map a
    # user-style chrom string (e.g. "1") to whatever this .hic file uses
    # (e.g. "chr1") when summing cis.
    chrom_norm_per_sample <- vector("list", length(hic_paths))
    hic_chroms_per_sample <- vector("list", length(hic_paths))
    for (i in seq_along(hic_paths)) {
        chroms_i <- strawr::readHicChroms(hic_paths[i])$name
        bp_res <- strawr::readHicBpResolutions(hic_paths[i])
        norms <- strawr::readHicNormTypes(hic_paths[i])
        if (!(resolution %in% bp_res)) {
            .stopNoCall(
                "Resolution ", resolution, " not available in ", hic_paths[i],
                " (available: ", paste(bp_res, collapse = ", "), ")"
            )
        }
        if (!("NONE" %in% norms)) {
            .stopNoCall("Normalization 'NONE' not available in ", hic_paths[i])
        }
        chrom_norm_per_sample[[i]] <- .make_chrom_normalizer(chroms_i)
        hic_chroms_per_sample[[i]] <- chroms_i
    }

    # Read and partition loops
    per_sample_loops <- lapply(loop_paths, .read_loops_bedpe)
    n_trans_per_sample <- vapply(
        per_sample_loops,
        function(d) sum(d$chrA != d$chrB),
        integer(1L)
    )
    per_sample_loops <- lapply(
        per_sample_loops,
        function(d) d[d$chrA == d$chrB, , drop = FALSE]
    )

    # Exact-coordinate consensus  # Canonicalise (smaller anchor first), build a string key per loop, then
    # count distinct samples per key. Keep keys with support >= min_samples.
    # No clustering, no jitter tolerance: two loops are the same only if
    # (chrA, startA, endA, chrB, startB, endB) match exactly.
    canonicalise <- function(d) {
        if (nrow(d) == 0L) {
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
    names(per_sample_loops) <- sample_ids

    per_sample_n <- vapply(per_sample_loops, nrow, integer(1L))
    if (sum(per_sample_n > 0L) < 2L) {
        .stopNoCall("Need at least 2 samples with at least one intra-chrom loop call")
    }
    if (any(per_sample_n == 0L)) {
        message(
            "Dropping ", sum(per_sample_n == 0L),
            " sample(s) with no usable loops: ",
            paste(sample_ids[per_sample_n == 0L], collapse = ", ")
        )
    }

    message(
        "[", format(Sys.time()), "] exact-coord consensus across ",
        sum(per_sample_n > 0L), " samples ..."
    )
    # Build a long data.frame: one row per per-sample call, with key + sample
    call_chunks <- lapply(seq_along(per_sample_loops), function(si) {
        d <- per_sample_loops[[si]]
        if (nrow(d) == 0L) {
            return(NULL)
        }
        data.frame(
            key = paste(d$chrA, d$startA, d$endA, d$startB, d$endB, sep = ":"),
            chrA = d$chrA,
            startA = d$startA, endA = d$endA,
            startB = d$startB, endB = d$endB,
            sample = sample_ids[si],
            stringsAsFactors = FALSE
        )
    })
    all_calls <- do.call(rbind, call_chunks)
    # Distinct-sample support per key (count unique samples after de-duping
    # within-sample repeats)
    uniq_call <- unique(all_calls[, c("key", "sample")])
    support <- table(uniq_call$key)
    keep_keys <- names(support)[support >= min_samples]
    message(
        "[", format(Sys.time()), "] ",
        length(support), " unique loop tuples; ",
        length(keep_keys), " survive min_samples = ", min_samples
    )
    if (length(keep_keys) == 0L) {
        .stopNoCall("No consensus loops survived min_samples = ", min_samples)
    }

    # Pick the first occurrence of each kept key for its coords
    first_idx <- match(keep_keys, all_calls$key)
    cons <- all_calls[first_idx, c("chrA", "startA", "endA", "startB", "endB"),
        drop = FALSE
    ]
    # Stable order: chrom, anchorA start, anchorB start
    ord <- order(cons$chrA, cons$startA, cons$startB)
    cons <- cons[ord, , drop = FALSE]
    rownames(cons) <- NULL
    n_support <- as.integer(support[paste(cons$chrA, cons$startA, cons$endA,
        cons$startB, cons$endB,
        sep = ":"
    )])

    # Build GInteractions (1-based inclusive ranges)
    gi <- InteractionSet::GInteractions(
        GenomicRanges::GRanges(
            cons$chrA,
            IRanges::IRanges(cons$startA + 1L, cons$endA)
        ),
        GenomicRanges::GRanges(
            cons$chrA,
            IRanges::IRanges(cons$startB + 1L, cons$endB)
        )
    )

    # Stable loop_id from anchor coordinates (0-based BEDPE-style starts)
    loop_id <- sprintf(
        "L%07d_%s_%d_%d",
        seq_along(gi),
        cons$chrA,
        cons$startA,
        cons$startB
    )
    S4Vectors::mcols(gi)$loop_id <- loop_id
    S4Vectors::mcols(gi)$n_samples <- n_support

    rep_A <- InteractionSet::anchors(gi, type = "first")

    # cis: user-supplied if provided, else auto-computed from .hic
    cis_chroms <- unique(as.character(GenomicRanges::seqnames(rep_A))) # chroms with consensus loops; defines the
    # set of chroms cis is summed over (consistent
    # with what normalize will process)

    cis_vec <- if (!is.null(user_cis)) {
        if (anyNA(user_cis) || any(user_cis <= 0)) {
            .stopNoCall("`samples$cis` must be positive numeric with no NAs")
        }
        setNames(as.numeric(user_cis), sample_ids)
    } else {
        message(
            "[", format(Sys.time()), "] computing cis depth from .hic files ",
            "(no `cis` column in samplesheet; ", length(cis_chroms),
            " chroms x ", length(sample_ids), " samples) ..."
        )
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
        cis        = cis_vec,
        vcf_path   = vcf,
        resolution = resolution,
        phenotype  = matrix(numeric(0), 0, 0),
        meta       = meta
    )
    methods::validObject(obj)
    obj
}

# Internal: compute per-sample cis depth by summing all straw counts across
# the given chroms. Returns a named numeric vector (one entry per sample).
# Parallel unit of work: one (sample, chrom) straw pull + sum.
.compute_cis_from_hic <- function(hic_paths, sample_ids, chroms,
                                  resolution, chrom_norms,
                                  BPPARAM = BiocParallel::SerialParam()) {
    grid <- expand.grid(
        sample_i = seq_along(sample_ids),
        chrom = chroms,
        stringsAsFactors = FALSE,
        KEEP.OUT.ATTRS = FALSE
    )

    one <- function(k) {
        si <- grid$sample_i[k]
        chr <- grid$chrom[k]
        chr_in_hic <- chrom_norms[[si]](chr)
        dat <- tryCatch(
            strawr::straw(
                "NONE", hic_paths[si], chr_in_hic, chr_in_hic,
                "BP", resolution
            ),
            error = function(e) NULL
        )
        if (is.null(dat) || nrow(dat) == 0L) {
            return(0)
        }
        sum(as.numeric(dat$counts))
    }

    results <- BiocParallel::bplapply(seq_len(nrow(grid)), one, BPPARAM = BPPARAM)
    cis_vec <- setNames(numeric(length(sample_ids)), sample_ids)
    for (k in seq_len(nrow(grid))) {
        cis_vec[grid$sample_i[k]] <- cis_vec[grid$sample_i[k]] + results[[k]]
    }
    cis_vec
}
