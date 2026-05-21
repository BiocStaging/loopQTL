# Light end-to-end check: build a 2-sample setup against strawr's toy .hic.
# Skips if strawr toy data isn't available in this install.

test_that("loopQTL_setup runs end-to-end against the strawr toy .hic", {
  hic <- toy_hic_path()
  td <- tempfile("loopQTL_toy_"); dir.create(td)

  bedpe <- write_toy_loops(td)

  # we don't have a tabix VCF on hand; build a minimal stub VCF so header
  # validation passes, then skip VCF-dependent assertions
  vcf_path <- file.path(td, "stub.vcf")
  cat("##fileformat=VCFv4.2\n",
      "##contig=<ID=1>\n",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
      "1\t1500000\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n",
      sep = "", file = vcf_path)
  vcf_gz <- paste0(vcf_path, ".gz")
  ok <- tryCatch({
    Rsamtools::bgzip(vcf_path, vcf_gz, overwrite = TRUE)
    Rsamtools::indexTabix(vcf_gz, format = "vcf")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) skip("Rsamtools bgzip/tabix not available in this environment")

  samples_df <- data.frame(
    Sample = c("A", "B"),
    hic    = c(hic, hic),
    loops  = c(bedpe, bedpe),
    stringsAsFactors = FALSE
  )

  obj <- loopQTL_setup(
    samples     = samples_df,
    vcf         = vcf_gz,
    resolution  = 250000L,
    min_samples = 1L
  )
  expect_s4_class(obj, "loopQTL")
  expect_equal(obj@samples, c("A", "B"))
  expect_equal(ncol(obj@contacts), 2L)
  expect_true(nrow(obj@contacts) >= 1L)
  expect_true(all(!is.na(obj@contacts)))   # NA -> 0 invariant
  expect_equal(names(obj@cis), c("A", "B"))
  expect_equal(obj@meta$cis_source, "computed_from_hic")
})

test_that("user-provided cis column overrides computed cis", {
  hic <- toy_hic_path()
  td <- tempfile("loopQTL_toy_"); dir.create(td)
  bedpe <- write_toy_loops(td)

  vcf_path <- file.path(td, "stub.vcf")
  cat("##fileformat=VCFv4.2\n",
      "##contig=<ID=1>\n",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
      "1\t1500000\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n",
      sep = "", file = vcf_path)
  vcf_gz <- paste0(vcf_path, ".gz")
  ok <- tryCatch({
    Rsamtools::bgzip(vcf_path, vcf_gz, overwrite = TRUE)
    Rsamtools::indexTabix(vcf_gz, format = "vcf")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) skip("Rsamtools bgzip/tabix not available in this environment")

  samples_df <- data.frame(
    Sample = c("A", "B"),
    hic    = c(hic, hic),
    loops  = c(bedpe, bedpe),
    cis    = c(1e7, 2e7),
    stringsAsFactors = FALSE
  )
  obj <- loopQTL_setup(samples_df, vcf_gz, resolution = 250000L, min_samples = 1L)
  expect_equal(unname(obj@cis), c(1e7, 2e7))
  expect_equal(obj@meta$cis_source, "user_provided")
})
