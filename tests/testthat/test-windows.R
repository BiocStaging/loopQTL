## Unit tests for the saveWindows feature.
#
## We avoid a real .hic by exercising the .compute_snp_windows helper directly
## with synthetic anchors and a small tabix-indexed VCF. The full end-to-end
## pass through loopQTLnormalize() is covered implicitly by the integration
## checks in test-setup.R / test-normalize.R when strawr's toy data is present.

make_tiny_vcf <- function(dir, snps) {
    ## snps: data.frame(chrom, pos, ref, alt, id)
    vcfPath <- file.path(dir, "tiny.vcf")
    con <- file(vcfPath, "w")
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
    vcfGzPath <- paste0(vcfPath, ".gz")
    Rsamtools::bgzip(vcfPath, vcfGzPath, overwrite = TRUE)
    Rsamtools::indexTabix(vcfGzPath, format = "vcf")
    vcfGzPath
}

## Two consensus loops on chr1 at 10 kb resolution:
##   loop 1: anchor A = bin 1_000_000, anchor B = bin 1_500_000
##   loop 2: anchor A = bin 1_200_000, anchor B = bin 2_000_000
## SNP s1 at 1_005_000  -> in anchor A of loop 1 only
## SNP s2 at 1_205_000  -> in anchor A of loop 2 only
## SNP s3 at 1_500_500  -> in anchor B of loop 1 only
## SNP s4 at 5_000_000  -> not in any anchor (should be dropped)
makeSynth <- function() {
    resolution <- 10000
    anchorA <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(1000001, 1200001),
            end = c(1010000, 1210000)
        )
    )
    anchorB <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(1500001, 2000001),
            end = c(1510000, 2010000)
        )
    )
    list(anchorA = anchorA, anchorB = anchorB, resolution = resolution)
}

test_that(".compute_snp_windows drops SNPs not in any anchor", {
    s <- makeSynth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcfGzPath <- make_tiny_vcf(td, data.frame(
        chrom = c("chr1", "chr1", "chr1", "chr1"),
        pos = c(1005000, 1205000, 1500500, 5000000),
        id = c("s1", "s2", "s3", "s4"),
        ref = c("A", "C", "G", "T"),
        alt = c("G", "T", "A", "C"),
        stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcfPath   = vcfGzPath,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snpWindow = 0,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_false(is.null(res))
    expect_equal(sort(res$snpsDF$snpID), c("s1", "s2", "s3"))
})

test_that(".compute_snp_windows returns NULL when no SNPs overlap any anchor", {
    s <- makeSynth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcfGzPath <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 5000000, id = "x",
        ref = "A", alt = "G", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcfPath   = vcfGzPath,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snpWindow = 0,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_null(res)
})

test_that("window extent is max(window, bounding span); cells form upper triangle", {
    s <- makeSynth()
    td <- tempfile("lqwin_")
    dir.create(td)
    ## SNP at the far end of loop-2's A anchor; partner B at 2_000_000.
    ## bounding span = 2_010_000 - 1_205_000 = 805_000 -> < window default 1e6
    ## so eff = 1e6, centered on midpoint ((1205000 + 2010000)/2) ~= 1_607_500
    vcfGzPath <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 1205000, id = "s2",
        ref = "C", alt = "T", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcfPath   = vcfGzPath,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snpWindow = 0,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_equal(nrow(res$snpsDF), 1)
    binLo <- res$snpsDF$binLo
    binHi <- res$snpsDF$binHi
    span <- binHi - binLo
    expect_gte(span, 1e6 - 2 * s$resolution) # rounded to bin grid
    ## Cell set is upper-triangle of binLo..binHi
    cc <- res$cells_by_chrom[["chr1"]]
    expect_true(all(cc$binI <= cc$binJ))
    expect_true(all(cc$binI >= binLo & cc$binJ <= binHi))
    ng <- (binHi - binLo) / s$resolution + 1
    expect_equal(nrow(cc), as.integer(ng * (ng + 1) / 2))
})

test_that("when bounding span exceeds window, window expands to cover it", {
    ## Build a synthetic far-apart loop pair: anchor A near 1Mb, anchor B near 5Mb
    ## span ~= 4Mb, > 1Mb default window -> eff should equal span (no shrinkage)
    resolution <- 10000
    anchorA <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 1000001, end = 1010000)
    )
    anchorB <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges   = IRanges::IRanges(start = 5000001, end = 5010000)
    )
    td <- tempfile("lqwin_")
    dir.create(td)
    vcfGzPath <- make_tiny_vcf(td, data.frame(
        chrom = "chr1", pos = 1005000, id = "s",
        ref = "A", alt = "G", stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcfPath   = vcfGzPath,
        anchorA    = anchorA,
        anchorB    = anchorB,
        snpWindow = 0,
        window     = 1e6, # smaller than bounding span ~4Mb
        resolution = resolution
    )
    expect_equal(nrow(res$snpsDF), 1)
    span <- res$snpsDF$binHi - res$snpsDF$binLo
    expect_gte(span, 4e6 - 2 * resolution) # at least the bounding extent
})

test_that("SNPs in different anchors of the same loop produce overlapping cells (deduplicated)", {
    s <- makeSynth()
    td <- tempfile("lqwin_")
    dir.create(td)
    vcfGzPath <- make_tiny_vcf(td, data.frame(
        chrom = c("chr1", "chr1"),
        pos = c(1005000, 1500500), # SNP s1 in anchor A; s3 in anchor B of loop 1
        id = c("s1", "s3"),
        ref = c("A", "G"),
        alt = c("G", "A"),
        stringsAsFactors = FALSE
    ))
    res <- loopQTL:::.compute_snp_windows(
        vcfPath   = vcfGzPath,
        anchorA    = s$anchorA,
        anchorB    = s$anchorB,
        snpWindow = 0,
        window     = 1e6,
        resolution = s$resolution
    )
    expect_equal(nrow(res$snpsDF), 2)
    ## Both SNPs live on chr1. Their union of cells is dedup'd.
    cc <- res$cells_by_chrom[["chr1"]]
    expect_true(!anyDuplicated(paste(cc$binI, cc$binJ, sep = "_")))
    ## Each SNP's local indices should point into cc, no NAs
    expect_false(anyNA(res$snp_cells_local[[1]]))
    expect_false(anyNA(res$snp_cells_local[[2]]))
})
