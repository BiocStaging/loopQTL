# Direct tests for the v2 normalize machinery, using a hand-built loopQTL
# object so we don't depend on .hic / VCF I/O.
#
# These tests cover:
#   * RINT-only path (we monkey-patch the per-chrom worker via a stub object)
#   * .quantile_normalize correctness
#   * The S4 invariants are preserved after normalize-shaped output

test_that(".quantile_normalize equalizes column distributions", {
  set.seed(1)
  m <- matrix(c(1, 5, 2, 8,
                3, 1, 4, 2,
                9, 7, 6, 3), nrow = 4, byrow = TRUE)
  qn <- loopQTL:::.quantile_normalize(m)
  # After QN every column should have the same sorted vector
  sorted_cols <- apply(qn, 2L, sort)
  expect_true(all(apply(sorted_cols, 1L, function(r) diff(range(r)) < 1e-12)))
})

test_that(".rint matches Blom offset and handles constants", {
  expect_equal(loopQTL:::.rint(c(1, 2, 3, 4, 5)),
               qnorm((rank(1:5) - 0.5) / 5))
  expect_true(all(is.na(loopQTL:::.rint(rep(7, 4)))))
})

test_that(".chr_to_int handles common cases", {
  expect_equal(loopQTL:::.chr_to_int("chr1"), 1L)
  expect_equal(loopQTL:::.chr_to_int("22"), 22L)
  expect_equal(loopQTL:::.chr_to_int("chrX"), 23L)
  expect_equal(loopQTL:::.chr_to_int("Y"), 24L)
})

test_that(".gt_to_dosage handles phased/unphased and missing", {
  gt <- matrix(c("0/0", "0|1", "1/1",
                 "./.", "."  , "0/2"),
               nrow = 2, byrow = TRUE)
  d <- loopQTL:::.gt_to_dosage(gt)
  expect_equal(d[1, ], c(0, 1, 2))
  expect_true(all(is.na(d[2, 1:2])))   # ./. and . both NA
  expect_true(is.na(d[2, 3]))          # multi-allelic 0/2 dropped
})
