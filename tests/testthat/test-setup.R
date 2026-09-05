## End-to-end smoke test of loopQTLsetup against strawr's toy .hic.
## Skips if strawr toy data isn't available in this install.
#
## v3 semantics:
##   - setup does NOT pull straw contacts; obj@contacts is empty by default
##   - obj@cis is empty unless samplesheet provides a `cis` column
##   - obj@metadata$hicPaths is populated (loopQTLnormalize needs it)

test_that("loopQTLsetup builds consensus loops without pulling contacts", {
    hic <- toy_hic_path()
    td <- tempfile("loopQTL_toy_")
    dir.create(td)
    bedpe <- write_toy_loops(td)

    vcfPath <- file.path(td, "stub.vcf")
    cat("##fileformat=VCFv4.2\n",
        "##contig=<ID=1>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        "1\t1500000\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n",
        sep = "", file = vcfPath
    )
    vcfGzPath <- paste0(vcfPath, ".gz")
    ok <- tryCatch(
        {
            Rsamtools::bgzip(vcfPath, vcfGzPath, overwrite = TRUE)
            Rsamtools::indexTabix(vcfGzPath, format = "vcf")
            TRUE
        },
        error = function(e) FALSE
    )
    if (!ok) skip("Rsamtools bgzip/tabix not available in this environment")

    samplesDF <- data.frame(
        Sample = c("A", "B"),
        hic = c(hic, hic),
        loops = c(bedpe, bedpe),
        stringsAsFactors = FALSE
    )

    obj <- loopQTLsetup(
        samples = samplesDF,
        vcf = vcfGzPath,
        resolution = 250000,
        minSamples = 1
    )
    expect_s4_class(obj, "loopQTL")
    expect_equal(obj@samples, c("A", "B"))

    ## Slim setup: contacts empty, cis empty (no samplesheet cis column)
    expect_equal(length(obj@contacts), 0)
    expect_equal(length(obj@cis), 0)
    expect_equal(obj@metadata$cisSource, "notProvided")

    ## hicPaths recorded for normalize
    expect_equal(unname(obj@metadata$hicPaths), c(hic, hic))
    expect_equal(names(obj@metadata$hicPaths), c("A", "B"))

    ## Consensus loops built
    expect_true(length(obj@loops) >= 1)
    expect_true(all(grepl("^L\\d{7}_", S4Vectors::mcols(obj@loops)$loopID)))
})

test_that("samplesheet cis column populates obj@cis", {
    hic <- toy_hic_path()
    td <- tempfile("loopQTL_toy_")
    dir.create(td)
    bedpe <- write_toy_loops(td)

    vcfPath <- file.path(td, "stub.vcf")
    cat("##fileformat=VCFv4.2\n",
        "##contig=<ID=1>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        "1\t1500000\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n",
        sep = "", file = vcfPath
    )
    vcfGzPath <- paste0(vcfPath, ".gz")
    ok <- tryCatch(
        {
            Rsamtools::bgzip(vcfPath, vcfGzPath, overwrite = TRUE)
            Rsamtools::indexTabix(vcfGzPath, format = "vcf")
            TRUE
        },
        error = function(e) FALSE
    )
    if (!ok) skip("Rsamtools bgzip/tabix not available in this environment")

    samplesDF <- data.frame(
        Sample = c("A", "B"),
        hic = c(hic, hic),
        loops = c(bedpe, bedpe),
        cis = c(1e7, 2e7),
        stringsAsFactors = FALSE
    )
    obj <- loopQTLsetup(samplesDF, vcfGzPath, resolution = 250000, minSamples = 1)
    expect_equal(unname(obj@cis), c(1e7, 2e7))
    expect_equal(names(obj@cis), c("A", "B"))
    expect_equal(obj@metadata$cisSource, "userProvided")
})
