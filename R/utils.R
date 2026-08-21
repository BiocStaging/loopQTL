## Internal function to display a warning without showing the function call
.warningNoCall <- function(...) {
    warning(..., call. = FALSE)
}

## Internal function to stop execution without showing the function call
.stopNoCall <- function(...) {
    stop(..., call. = FALSE)
}

## Internal function to display an error like stop() would, but without stopping
## execution or showing the function call
.displayError <- function(...) {
    try(stop(..., call. = FALSE))
}

# Bin a 1-based BEDPE-style anchor (start, end) to its analysis-bin start (0-based),
# matching strawr's binning. straw returns bin starts as the lower coordinate of
# the bin in 0-based half-open BP space, so bin = floor(midpoint / resolution) * resolution.
.bin_anchor <- function(start, end, resolution) {
    mid <- (as.numeric(start) + as.numeric(end)) / 2
    as.integer(floor(mid / resolution) * resolution)
}

# Detect strawr's chromosome naming style ("1" vs "chr1") from the .hic header.
# Returns a function that maps a user-supplied chrom string (either style) to the
# style this .hic expects.
.make_chrom_normalizer <- function(hic_chroms) {
    uses_chr_prefix <- any(grepl("^chr", hic_chroms))
    function(x) {
        if (uses_chr_prefix) {
            ifelse(grepl("^chr", x), x, paste0("chr", x))
        } else {
            sub("^chr", "", x)
        }
    }
}

# RINT with Blom offset, NA-aware. Operates on a single vector (one loop's
# values across samples).
.rint <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    n <- sum(!is.na(r))
    if (n < 2L) {
        return(rep(NA_real_, length(x)))
    }
    qnorm((r - 0.5) / n)
}

# Apply .rint() row-wise to a (loops x samples) matrix. Rows whose values are
# all-NA, fewer than 2 non-NA, or have zero variance are returned as all-NA
# (RINT on a constant produces NaN/Inf, which would poison downstream
# regressions). dimnames are preserved.
.rint_matrix <- function(mat) {
    if (length(mat) == 0L) {
        return(mat)
    }
    out <- matrix(NA_real_,
        nrow = nrow(mat), ncol = ncol(mat),
        dimnames = dimnames(mat)
    )
    row_sd <- apply(mat, 1L, stats::sd, na.rm = TRUE)
    ok <- is.finite(row_sd) & row_sd > 0
    if (any(ok)) {
        good <- which(ok)
        out[good, ] <- t(apply(mat[good, , drop = FALSE], 1L, .rint))
    }
    out
}

# Internal: compute per-SNP window cells for save_windows in loopQTLnormalize().
#
# For each VCF SNP that overlaps any consensus loop anchor (within `snp_window`
# bp flank), build the SNP's "window box" - a square (bin_lo..bin_hi)^2 region
# centered on the bounding interval of the SNP position and every partner
# anchor of every loop the SNP participates in, expanded symmetrically so its
# span is at least `window` bp. Returns the per-SNP metadata plus per-chrom
# deduplicated cell-coordinate tables that the BNBC per-chrom worker will
# index after .normalize_BNBC() emits its wide table.
#
# Returns:
#   list(
#     snps_df = data.frame(snp_id, chrom_orig, chrom_vcf, pos, ref, alt,
#                          bin_lo, bin_hi, n_loops_anchored),
#     cells_by_chrom = named list per chrom; each is data.frame(bin_i, bin_j)
#                      of upper-triangle cells (bp), deduplicated across SNPs
#                      on that chrom,
#     snp_cells_local = list (length nrow(snps_df)) of integer vectors of
#                      within-chrom row indices into cells_by_chrom[[chrom]]
#   )
# Returns NULL if no eligible SNPs are found (caller short-circuits).
.compute_snp_windows <- function(vcf_path, anchorA, anchorB,
                                 snp_window, window, resolution) {
    if (!file.exists(vcf_path)) {
        .stopNoCall("VCF not found: ", vcf_path)
    }

    # VCF chrom-naming normalizer
    tab <- Rsamtools::TabixFile(vcf_path)
    vcf_chroms <- Rsamtools::headerTabix(tab)$seqnames
    if (length(vcf_chroms) == 0L) {
        .stopNoCall("Could not read chromosome names from VCF tabix index")
    }
    vcf_norm <- .make_chrom_normalizer(vcf_chroms)

    # Build SNP target ranges = anchors padded by snp_window
    pad <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = vcf_norm(as.character(GenomicRanges::seqnames(gr))),
            ranges = IRanges::IRanges(
                start = pmax(1L, GenomicRanges::start(gr) - snp_window),
                end   = GenomicRanges::end(gr) + snp_window
            )
        )
    }
    targets <- GenomicRanges::reduce(c(pad(anchorA), pad(anchorB)))
    targets <- targets[as.character(GenomicRanges::seqnames(targets)) %in% vcf_chroms]
    if (length(targets) == 0L) {
        return(NULL)
    }

    message(
        "save_windows: tabix-reading SNPs at ",
        length(targets), " anchor regions (snp_window = ", snp_window, " bp) ..."
    )
    param <- VariantAnnotation::ScanVcfParam(
        which = targets,
        geno = NA_character_,
        info = NA_character_
    )
    vcf <- VariantAnnotation::readVcf(vcf_path, genome = "unknown", param = param)
    if (length(vcf) == 0L) {
        warning(
            "save_windows: no SNPs in any anchor +/- ",
            snp_window, " bp window"
        )
        return(NULL)
    }
    rr <- SummarizedExperiment::rowRanges(vcf)
    snp_chrom_vcf <- as.character(GenomicRanges::seqnames(rr))
    snp_pos <- as.integer(GenomicRanges::start(rr))
    snp_ref <- as.character(VariantAnnotation::ref(vcf))
    alt_list <- VariantAnnotation::alt(vcf)
    snp_alt <- vapply(
        alt_list,
        function(a) paste(as.character(a), collapse = ","),
        character(1L)
    )
    snp_id <- rownames(vcf)
    if (is.null(snp_id) || any(snp_id == "" | is.na(snp_id))) {
        snp_id <- sprintf("%s:%d_%s/%s", snp_chrom_vcf, snp_pos, snp_ref, snp_alt)
    }

    # Map SNPs to consensus loops (within snp_window flank)  # Convert anchors to vcf-style seqnames for overlap with the SNP GRanges.
    to_vcf_gr <- function(gr) {
        GenomicRanges::GRanges(
            seqnames = vcf_norm(as.character(GenomicRanges::seqnames(gr))),
            ranges = IRanges::IRanges(
                start = pmax(1L, GenomicRanges::start(gr) - snp_window),
                end   = GenomicRanges::end(gr) + snp_window
            )
        )
    }
    padded_A_vcf <- to_vcf_gr(anchorA)
    padded_B_vcf <- to_vcf_gr(anchorB)

    snp_gr <- GenomicRanges::GRanges(
        seqnames = snp_chrom_vcf,
        ranges   = IRanges::IRanges(start = snp_pos, end = snp_pos)
    )

    hitsA <- GenomicRanges::findOverlaps(snp_gr, padded_A_vcf)
    hitsB <- GenomicRanges::findOverlaps(snp_gr, padded_B_vcf)
    in_loops_by_snp <- vector("list", length(snp_gr))
    if (length(hitsA) > 0L) {
        spl <- split(S4Vectors::subjectHits(hitsA), S4Vectors::queryHits(hitsA))
        for (s in names(spl)) {
            in_loops_by_snp[[as.integer(s)]] <- c(in_loops_by_snp[[as.integer(s)]], spl[[s]])
        }
    }
    if (length(hitsB) > 0L) {
        spl <- split(S4Vectors::subjectHits(hitsB), S4Vectors::queryHits(hitsB))
        for (s in names(spl)) {
            in_loops_by_snp[[as.integer(s)]] <- c(in_loops_by_snp[[as.integer(s)]], spl[[s]])
        }
    }
    keep_snp <- vapply(in_loops_by_snp, function(x) length(x) > 0L, logical(1L))
    if (!any(keep_snp)) {
        .warningNoCall(
            "save_windows: no SNPs overlap any consensus loop anchor"
        )
        return(NULL)
    }
    snp_idx_keep <- which(keep_snp)
    in_loops_by_snp <- in_loops_by_snp[snp_idx_keep]

    snp_id <- snp_id[snp_idx_keep]
    snp_chrom_vcf <- snp_chrom_vcf[snp_idx_keep]
    snp_pos <- snp_pos[snp_idx_keep]
    snp_ref <- snp_ref[snp_idx_keep]
    snp_alt <- snp_alt[snp_idx_keep]

    # Native (anchor-style) chrom for output. Inverse of vcf_norm is sample-specific;
    # use anchor seqnames at the SNP's chrom directly.
    anchor_chrom_native <- as.character(GenomicRanges::seqnames(anchorA))
    anchor_chrom_vcf <- vcf_norm(anchor_chrom_native)
    snp_chrom_orig <- vapply(snp_chrom_vcf, function(c_vcf) {
        m <- which(anchor_chrom_vcf == c_vcf)[1L]
        if (is.na(m)) c_vcf else anchor_chrom_native[m]
    }, character(1L))

    # Per-SNP bounding interval -> bin_lo/bin_hi
    anchorA_start_native <- GenomicRanges::start(anchorA)
    anchorA_end_native <- GenomicRanges::end(anchorA)
    anchorB_start_native <- GenomicRanges::start(anchorB)
    anchorB_end_native <- GenomicRanges::end(anchorB)

    n_snp <- length(snp_id)
    bin_lo <- integer(n_snp)
    bin_hi <- integer(n_snp)
    n_loops <- integer(n_snp)

    for (i in seq_len(n_snp)) {
        li <- unique(in_loops_by_snp[[i]])
        n_loops[i] <- length(li)
        # Bounding interval: SNP pos U all involved anchors (A and B)
        starts <- c(snp_pos[i], anchorA_start_native[li], anchorB_start_native[li])
        ends <- c(snp_pos[i], anchorA_end_native[li], anchorB_end_native[li])
        b_lo <- min(starts)
        b_hi <- max(ends)
        span <- b_hi - b_lo
        if (span < window) {
            center <- (b_lo + b_hi) / 2
            half <- window / 2
            b_lo <- center - half
            b_hi <- center + half
        }
        bin_lo[i] <- as.integer(floor(b_lo / resolution) * resolution)
        bin_hi[i] <- as.integer(floor(b_hi / resolution) * resolution)
        if (bin_lo[i] < 0L) bin_lo[i] <- 0L
    }

    snps_df <- data.frame(
        snp_id = snp_id,
        chrom_orig = snp_chrom_orig,
        chrom_vcf = snp_chrom_vcf,
        pos = snp_pos,
        ref = snp_ref,
        alt = snp_alt,
        bin_lo = bin_lo,
        bin_hi = bin_hi,
        n_loops_anchored = n_loops,
        stringsAsFactors = FALSE
    )

    # Per-chrom dedup of cells  # cell = (bin_i, bin_j) with bin_lo <= bin_i <= bin_j <= bin_hi (upper tri).
    # Group SNPs by their native (anchor-style) chrom, since BNBC emits in that style.
    cells_by_chrom <- list()
    snp_cells_local <- vector("list", n_snp)

    chrom_groups <- split(seq_len(n_snp), snps_df$chrom_orig)
    for (chr in names(chrom_groups)) {
        idxs <- chrom_groups[[chr]]
        # Build a master set of cells, dedupe via paste-key
        cell_chunks <- vector("list", length(idxs))
        for (k in seq_along(idxs)) {
            s <- idxs[k]
            bins <- seq.int(snps_df$bin_lo[s], snps_df$bin_hi[s], by = resolution)
            ng <- length(bins)
            if (ng == 0L) {
                cell_chunks[[k]] <- data.frame(
                    bin_i = integer(0),
                    bin_j = integer(0)
                )
                next
            }
            # All (bin_i, bin_j) with bin_i <= bin_j, upper triangle inclusive
            ii <- rep(bins, times = ng:1)
            jj <- unlist(lapply(seq_len(ng), function(z) bins[z:ng]))
            cell_chunks[[k]] <- data.frame(
                bin_i = as.integer(ii),
                bin_j = as.integer(jj)
            )
        }
        combined <- do.call(rbind, cell_chunks)
        keys <- paste(combined$bin_i, combined$bin_j, sep = "_")
        uniq_idx <- !duplicated(keys)
        chrom_cells <- combined[uniq_idx, , drop = FALSE]
        rownames(chrom_cells) <- NULL
        chrom_cells_keys <- keys[uniq_idx]
        cells_by_chrom[[chr]] <- chrom_cells

        # For each SNP, look up its cells' local indices
        for (k in seq_along(idxs)) {
            s <- idxs[k]
            ch <- cell_chunks[[k]]
            if (nrow(ch) == 0L) {
                snp_cells_local[[s]] <- integer(0)
                next
            }
            ck <- paste(ch$bin_i, ch$bin_j, sep = "_")
            snp_cells_local[[s]] <- match(ck, chrom_cells_keys)
        }
    }

    list(
        snps_df         = snps_df,
        cells_by_chrom  = cells_by_chrom,
        snp_cells_local = snp_cells_local
    )
}

# Read a single BEDPE-like loop file. Tolerant of "#chr" headers and of
# files with or without a `loop_id` column.
.read_loops_bedpe <- function(path) {
    if (!file.exists(path)) {
        .stopNoCall("Loop file not found: ", path)
    }
    hdr <- readLines(path, n = 1L, warn = FALSE)
    has_header <- grepl("^#?chr", hdr, ignore.case = TRUE) &&
        grepl("\\bstart", hdr, ignore.case = TRUE)
    df <- utils::read.table(
        path,
        header = has_header, sep = "\t",
        comment.char = if (has_header) "" else "#",
        stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!has_header) {
        # assume positional: chrA startA endA chrB startB endB [loop_id]
        if (ncol(df) < 6L) {
            .stopNoCall("Loop file ", path, " must have at least 6 columns")
        }
        base_names <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
        if (ncol(df) == 6L) {
            names(df) <- base_names
        } else {
            names(df) <- c(base_names, paste0("V", seq_len(ncol(df) - 6L) + 6L))
        }
    } else {
        nm <- names(df)
        if (nm[1] %in% c("#chr1", "#chrA")) nm[1] <- "chrA"
        # canonicalize common variants
        nm <- sub("^chr1$", "chrA", nm)
        nm <- sub("^x1$", "startA", nm, ignore.case = TRUE)
        nm <- sub("^y1$", "endA", nm, ignore.case = TRUE)
        nm <- sub("^chr2$", "chrB", nm)
        nm <- sub("^x2$", "startB", nm, ignore.case = TRUE)
        nm <- sub("^y2$", "endB", nm, ignore.case = TRUE)
        names(df) <- nm
    }
    need <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
    miss <- setdiff(need, names(df))
    if (length(miss) > 0L) {
        .stopNoCall(
            "Loop file ", path, " missing required columns: ",
            paste(miss, collapse = ", ")
        )
    }
    df$chrA <- as.character(df$chrA)
    df$chrB <- as.character(df$chrB)
    df$startA <- as.integer(df$startA)
    df$endA <- as.integer(df$endA)
    df$startB <- as.integer(df$startB)
    df$endB <- as.integer(df$endB)
    df
}
