test_that(".binAnchor floors midpoint to bin start", {
    expect_equal(loopQTL:::.binAnchor(1, 250000, 250000), 125000)
    expect_equal(loopQTL:::.binAnchor(250000, 500000, 250000), 375000 %/% 250000 * 250000)
    expect_equal(loopQTL:::.binAnchor(0, 250000, 250000), 125000 %/% 250000 * 250000)
})

test_that(".rint matches Blom offset", {
    x <- c(1, 2, 3, 4, 5)
    r <- loopQTL:::.rint(x)
    expect_equal(r, qnorm((rank(x) - 0.5) / length(x)))
})

test_that(".rint returns NA on constant input", {
    expect_true(all(is.na(loopQTL:::.rint(rep(1, 5)))))
})

test_that("chrom normalizer handles both styles", {
    fChr <- loopQTL:::.make_chrom_normalizer(c("chr1", "chr2", "chrX"))
    expect_equal(fChr("1"), "chr1")
    expect_equal(fChr("chr1"), "chr1")
    fNochr <- loopQTL:::.make_chrom_normalizer(c("1", "2", "X"))
    expect_equal(fNochr("chr1"), "1")
    expect_equal(fNochr("1"), "1")
})
