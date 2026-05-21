# loopQTL

Anchor-restricted loop QTL mapping from Hi-C — a single-environment, pure-R package for testing whether genotype at a chromatin loop's anchors associates with the loop's contact strength across samples.

Unlike peak-centric, window-scan tools (e.g. QTLtools), `loopQTL` only tests SNPs that fall in the loop's **anchor bins**, because most window-scan SNPs are biologically meaningless for a given loop. Intra-chromosomal loops only.

## Pipeline

```
loopQTL_setup()    -> reads .hic per sample, builds bin-matched consensus loops,
                      builds (loops x samples) contact matrix
loopQTL_normalize() -> DESeq2 -> log2 -> residualize log10(cis+1) (+ user covariates)
                      -> RINT  (single RINT, owned here)
run_QTL()          -> targeted tabix VCF read of anchor SNPs, lm() per (loop, SNP)
                      -> tidy results data.frame
```

A single `loopQTL` S4 object is threaded through all three functions. The canonical sample order is an invariant; matrices stay aligned.

## Install

The package depends on a handful of Bioconductor packages. From R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "strawr", "DESeq2", "GenomicRanges", "IRanges", "InteractionSet",
  "S4Vectors", "VariantAnnotation", "Rsamtools",
  "SummarizedExperiment", "BiocParallel"
))

# install loopQTL from this repo:
# devtools::install_github("tim-wuuu/loopQTL")
# or, from a local clone:
# devtools::install_local("path/to/loopQTL")
```

For iterative development without installing:

```r
devtools::load_all("path/to/loopQTL")
devtools::test()
```

## Minimal example

```r
library(loopQTL)

samples_df <- data.frame(
  Sample = c("S1", "S2", "S3"),
  hic    = c("S1.hic", "S2.hic", "S3.hic"),
  loops  = c("S1_loops.bedpe", "S2_loops.bedpe", "S3_loops.bedpe"),
  # optional: precomputed cis-contact totals; otherwise computed from .hic
  cis    = c(8.1e7, 7.4e7, 9.2e7),
  stringsAsFactors = FALSE
)

obj <- loopQTL_setup(
  samples     = samples_df,
  vcf         = "cohort.vcf.gz",     # bgzipped + tabix-indexed
  resolution  = 10000L,              # bp; must exist in every .hic
  min_samples = 2L,                  # consensus filter
  BPPARAM     = BiocParallel::MulticoreParam(workers = 8)
)

obj <- loopQTL_normalize(obj)
# optional extra covariates (e.g. genotype PCs), samples x k:
# obj <- loopQTL_normalize(obj, covariates = pcs)

results <- run_QTL(obj, window = 0L)   # window = bp flank around anchor bin
head(results)
#   loop_id   snp_id          chrom  pos      ref alt anchor   beta   se     t      p   n_used  maf
```

## Design decisions

These are settled — see the [design spec](../loopQTL_package_spec.md) for the full rationale:

- **Anchor-restricted, not window-scan.** SNPs are only tested if they fall in (or within `window` of) the loop's anchor bins.
- **Intra-chromosomal loops only.** Trans loops are dropped at setup; count recorded in `meta`.
- **Bin-match consensus.** Two loops are "the same" if both anchors fall in the same `resolution` bin pair.
- **Uniform NA -> 0.** Any absent contact lookup becomes 0. There are no NAs in the contact matrix after setup.
- **Single RINT, fixed order.** DESeq2 -> log2 -> residualize -> RINT, applied exactly once, inside `loopQTL_normalize()`.
- **`run_QTL()` is bare.** No transformation, no covariate adjustment, no FDR, no lead-SNP selection. Returns the full per-(loop, SNP) table; significance calling is up to you.

## Parallelism

The unit of parallel work is one `(sample, chromosome)` extraction in `loopQTL_setup()`. Dispatched via `BiocParallel::bplapply` — the backend is your choice via `BPPARAM`:

```r
# workstation:
BPPARAM = BiocParallel::MulticoreParam(workers = 16)

# SLURM cluster:
BPPARAM = BiocParallel::BatchtoolsParam(workers = 200, cluster = "slurm",
                                        template = "slurm.tmpl")
```

## Testing

A small `testthat` suite uses `strawr`'s bundled `test.hic` plus an on-the-fly stub VCF for an end-to-end smoke test. Run with:

```r
devtools::test()
# or
R CMD check .
```

## License

MIT. See [LICENSE](LICENSE).
