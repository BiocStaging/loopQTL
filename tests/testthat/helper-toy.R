# Helpers for tests. The toy .hic comes from strawr; everything else is built on the fly.

toy_hic_path <- function() {
  p <- system.file("extdata", "test.hic", package = "strawr")
  if (!nzchar(p)) skip("strawr toy .hic not available")
  p
}

# Build a tiny BEDPE loop file derived from real entries in the toy .hic so
# extraction has something to find. Two loops on chr1 at 250kb resolution.
write_toy_loops <- function(dir, fname = "loops.bedpe", resolution = 250000L) {
  # straw uses bin starts; the toy file has chr "1" at 250000 bp resolution.
  df <- data.frame(
    chrA   = c("1", "1"),
    startA = c(1000000L, 5000000L),
    endA   = c(1000000L + resolution, 5000000L + resolution),
    chrB   = c("1", "1"),
    startB = c(2000000L, 6000000L),
    endB   = c(2000000L + resolution, 6000000L + resolution),
    stringsAsFactors = FALSE
  )
  path <- file.path(dir, fname)
  utils::write.table(df, path, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  path
}
