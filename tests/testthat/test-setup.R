# End-to-end smoke test of loopQTL_setup against strawr's toy .hic.
# Skips if strawr toy data isn't available in this install.
#
# v3 semantics:
#   - setup does NOT pull straw contacts; obj@contacts is empty by default
#   - obj@cis is empty unless samplesheet provides a `cis` column
#   - obj@meta$hic_paths is populated (loopQTL_normalize needs it)

test_that("loopQTL_setup builds consensus loops without pulling contacts", {
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

  # Slim setup: contacts empty, cis empty (no samplesheet cis column)
  expect_equal(length(obj@contacts), 0L)
  expect_equal(length(obj@cis), 0L)
  expect_equal(obj@meta$cis_source, "not_provided")

  # hic_paths recorded for normalize
  expect_equal(unname(obj@meta$hic_paths), c(hic, hic))
  expect_equal(names(obj@meta$hic_paths), c("A", "B"))

  # Consensus loops built
  expect_true(length(obj@loops) >= 1L)
  expect_true(all(grepl("^L\\d{7}_", S4Vectors::mcols(obj@loops)$loop_id)))
})

test_that("samplesheet cis column populates obj@cis", {
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
  expect_equal(names(obj@cis), c("A", "B"))
  expect_equal(obj@meta$cis_source, "user_provided")
})
