# Direct unit tests for the normalization math, constructed from a manual
# loopQTL object so we don't depend on .hic / VCF I/O.

make_fake_obj <- function(n_loops = 50, n_samples = 8, seed = 1) {
  set.seed(seed)
  samples_vec <- paste0("S", seq_len(n_samples))
  counts <- matrix(rpois(n_loops * n_samples, lambda = 30),
                   nrow = n_loops, ncol = n_samples,
                   dimnames = list(NULL, samples_vec))
  # induce a depth gradient
  depth <- seq(0.5, 1.5, length.out = n_samples)
  counts <- round(t(t(counts) * depth))
  storage.mode(counts) <- "double"
  cis_vec <- setNames(colSums(counts) * 10, samples_vec)

  gi <- InteractionSet::GInteractions(
    GenomicRanges::GRanges("1", IRanges::IRanges(start = seq_len(n_loops) * 1e5,
                                                  width = 1e4)),
    GenomicRanges::GRanges("1", IRanges::IRanges(start = seq_len(n_loops) * 1e5 + 5e5,
                                                  width = 1e4))
  )
  S4Vectors::mcols(gi)$loop_id <- sprintf("L%03d", seq_len(n_loops))

  new("loopQTL",
      samples    = samples_vec,
      loops      = gi,
      contacts   = counts,
      cis        = cis_vec,
      vcf_path   = tempfile(fileext = ".vcf.gz"),
      resolution = 10000L,
      phenotype  = matrix(numeric(0), 0, 0),
      meta       = list())
}

test_that("loopQTL_normalize fills phenotype with correct dims and names", {
  obj <- make_fake_obj()
  obj2 <- loopQTL_normalize(obj)
  expect_equal(dim(obj2@phenotype), dim(obj@contacts))
  expect_equal(colnames(obj2@phenotype), obj@samples)
})

test_that("RINT row-marginals are approximately standard normal", {
  obj <- loopQTL_normalize(make_fake_obj(n_loops = 200, n_samples = 20))
  # each row was RINT'd; its mean should be ~0 and sd ~1
  rm <- rowMeans(obj@phenotype, na.rm = TRUE)
  rs <- apply(obj@phenotype, 1, sd, na.rm = TRUE)
  expect_true(median(abs(rm)) < 1e-6)
  expect_true(median(rs, na.rm = TRUE) > 0.8)
})

test_that("zero-variance row becomes NA, count tallied in meta", {
  obj <- make_fake_obj()
  obj@contacts[1, ] <- 5   # constant row -> residuals constant -> NA after detection
  obj2 <- loopQTL_normalize(obj)
  expect_true(all(is.na(obj2@phenotype[1, ])))
  expect_true(obj2@meta$n_zero_var >= 1L)
})

test_that("extra covariate is accepted and applied", {
  obj <- make_fake_obj()
  extra <- data.frame(pc1 = rnorm(length(obj@samples)),
                      row.names = obj@samples)
  obj2 <- loopQTL_normalize(obj, covariates = extra)
  expect_true("pc1" %in% obj2@meta$covariates)
  expect_equal(dim(obj2@phenotype), dim(obj@contacts))
})

test_that("normalize order is fixed and recorded", {
  obj2 <- loopQTL_normalize(make_fake_obj())
  expect_equal(obj2@meta$normalize_order,
               c("DESeq2", "log2", "residualize", "RINT"))
})
