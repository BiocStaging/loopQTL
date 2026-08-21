## Unit tests for the save_windows feature.
#
## We avoid a real .hic by exercising the .compute_snp_windows helper directly
## with synthetic anchors and a small tabix-indexed VCF. The full end-to-end
## pass through loopQTLnormalize() is covered implicitly by the integration
## checks in test-setup.R / test-normalize.R when strawr's toy data is present.

make_tiny_vcf <- function(dir, snps) {
    ## snps: data.frame(chrom, pos, ref, alt, id)
    vcf_path <- file.path(dir, "tiny.vcf")
    con <- file(vcf_path, "w")
    on.exit(close(con), add = TRUE)
    cat("##fileformat=VCFv4.2\n", file = con)
    for (chr in unique(snps$chrom)) {
        cat("##contig=<ID=", chr, ">\n", sep = "", file = con)
    }
    cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n", file = con)
    for (i in seq_len(nrow(snps))) {
        cat(snps$chrom[i], "\t", snps$pos[i], "\t", snps$id[i], "\t",
            snps$ref[i], "\t", snps$alt[i],
            "\t.\tPASS\t.\tGT\t0/1\t1/1\n",
            sep = "", file = con
        )
    }
    close(con)
    on.exit()
    vcf_gz <- paste0(vcf_path, ".gz")
    Rsamtools::bgzip(vcf_path, vcf_gz, overwrite = TRUE)
    Rsamtools::indexTabix(vcf_gz, format = "vcf")
    vcf_gz
}

## Two consensus loops on chr1 at 10 kb resolution:
##   loop 1: anchor A = bin 1_000_000, anchor B = bin 1_500_000
##   loop 2: anchor A = bin 1_200_000, anchor B = bin 2_000_000
## SNP s1 at 1_005_000  -> in anchor A of loop 1 only
## SNP s2 at 1_205_000  -> in anchor A of loop 2 only
## SNP s3 at 1_500_500  -> in anchor B of loop 1 only
## SNP s4 at 5_000_000  -> not in any anchor (should be dropped)
make_synth <- function() {
    resolution <- 10000L
    anchorA <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(1000001L, 1200001L),
            end = c(1010000L, 1210000L)
        )
    )
    anchorB <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(1500001L, 2000001L),
            end = c(1510000L, 2010000L)
        )
    )
    list(anchorA = anchorA, anchorB = anchorB, resolution = resolution)
}

test_that(".compute_snp_windows drops SNPs not in any anchor", {
    s <- make_synth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcf_gz <- make_tiny_vcf(td, data.frame(
        chrom = c("chr1", "chr1", "chr1", "chr1"),
        pos = c(1005000L, 1205000L, 1500500L, 5000000L),
        id = c("s1", "s2", "s3", "s4"),
        ref = c("A", "C", "G", "T"),
        alt = c("G", "T", "A", "C"),
        stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcf_path   = vcf_gz,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snp_window = 0L,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_false(is.null(res))
    expect_equal(sort(res$snps_df$snp_id), c("s1", "s2", "s3"))
})

test_that(".compute_snp_windows returns NULL when no SNPs overlap any anchor", {
    s <- make_synth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcf_gz <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 5000000L, id = "x",
        ref = "A", alt = "G", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcf_path   = vcf_gz,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snp_window = 0L,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_null(res)
})

test_that("window extent is max(window, bounding span); cells form upper triangle", {
    s <- make_synth()
    td <- tempfile("lqwin_")
    dir.create(td)
    ## SNP at the far end of loop-2's A anchor; partner B at 2_000_000.
    ## bounding span = 2_010_000 - 1_205_000 = 805_000 -> < window default 1e6
    ## so eff = 1e6, centered on midpoint ((1205000 + 2010000)/2) ~= 1_607_500
    vcf_gz <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 1205000L, id = "s2",
        ref = "C", alt = "T", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcf_path   = vcf_gz,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snp_window = 0L,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_equal(nrow(res$snps_df), 1L)
    bin_lo <- res$snps_df$bin_lo
    bin_hi <- res$snps_df$bin_hi
    span <- bin_hi - bin_lo
    expect_gte(span, 1e6 - 2 * s$resolution) # rounded to bin grid
    ## Cell set is upper-triangle of bin_lo..bin_hi
    cc <- res$cells_by_chrom[["chr1"]]
    expect_true(all(cc$bin_i <= cc$bin_j))
    expect_true(all(cc$bin_i >= bin_lo & cc$bin_j <= bin_hi))
    ng <- (bin_hi - bin_lo) / s$resolution + 1L
    expect_equal(nrow(cc), as.integer(ng * (ng + 1L) / 2L))
})

test_that("when bounding span exceeds window, window expands to cover it", {
    ## Build a synthetic far-apart loop pair: anchor A near 1Mb, anchor B near 5Mb
    ## span ~= 4Mb, > 1Mb default window -> eff should equal span (no shrinkage)
    resolution <- 10000L
    anchorA <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 1000001L, end = 1010000L)
    )
    anchorB <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 5000001L, end = 5010000L)
    )
    td <- tempfile("lqwin_")
    dir.create(td)
    vcf_gz <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 1005000L, id = "s",
        ref = "A", alt = "G", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcf_path   = vcf_gz,
        anchorA    = anchorA,
        anchorB    = anchorB,
        snp_window = 0L,
        window     = 1e6, # smaller than bounding span ~4Mb
        resolution = resolution
    )
    expect_equal(nrow(res$snps_df), 1L)
    span <- res$snps_df$bin_hi - res$snps_df$bin_lo
    expect_gte(span, 4e6 - 2 * resolution) # at least the bounding extent
})

test_that("SNPs in different anchors of the same loop produce overlapping cells (deduplicated)", {
    s <- make_synth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcf_gz <- make_tiny_vcf(td, data.frame(
        chrom = c("chr1", "chr1"),
        pos = c(1005000L, 1500500L), # SNP s1 in anchor A; s3 in anchor B of loop 1
        id = c("s1", "s3"),
        ref = c("A", "G"),
        alt = c("G", "A"),
        stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcf_path   = vcf_gz,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snp_window = 0L,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_equal(nrow(res$snps_df), 2L)
    ## Both SNPs live on chr1. Their union of cells is dedup'd.
    cc <- res$cells_by_chrom[["chr1"]]
    expect_true(!anyDuplicated(paste(cc$bin_i, cc$bin_j, sep = "_")))
    ## Each SNP's local indices should point into cc, no NAs
    expect_false(anyNA(res$snp_cells_local[[1]]))
    expect_false(anyNA(res$snp_cells_local[[2]]))
})
