# loopQTL

An R package for testing whether genotype at a chromatin loop's anchors associates with the loop's contact strength across samples.

Unlike peak-centric, window-scan tools (e.g. QTLtools), `loopQTL` only tests SNPs that fall in the loop's **anchor bins**, because most window-scan SNPs are biologically meaningless for a given loop. Intra-chromosomal loops only. Meanwhile, it considers the total cis reads depth for the normalization.

## Pipeline

```
loopQTL_setup()    -> validates VCF + .hic headers, reads loop BEDPEs,
                      builds bin-matched consensus loops. No contact extraction
                      (fast; iterate on the samplesheet freely).
loopQTL_normalize() -> per-chrom whole-genome straw pull (parallel),
                      multiHiCcompare fastlo() or BNBC band-wise
                      quantile normalization, extract at consensus loops,
                      optional RINT (default TRUE)
loopQTL_choose_K()  -> PCAForQTL elbow/BE selection of K phenotype PCs (helper)
run_QTL()          -> targeted tabix VCF read of anchor SNPs, MatrixEQTL run
                      once per anchor, merge + label anchor (A/B/both)
                      -> tidy results data.frame
```

A single `loopQTL` S4 object is threaded through. The canonical sample order is an invariant; matrices stay aligned.

## Install

The package depends on a handful of Bioconductor packages. From R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "strawr", "multiHiCcompare", "preprocessCore",
  "GenomicRanges", "IRanges", "InteractionSet",
  "S4Vectors", "VariantAnnotation", "Rsamtools",
  "SummarizedExperiment", "BiocParallel"
))
install.packages(c("MatrixEQTL", "data.table"))
# PCAForQTL is GitHub-only:
remotes::install_github("heatherjzhou/PCAForQTL")

# install loopQTL from this repo:
# remotes::install_github("tim-wuuu/loopQTL")
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
  min_samples = 2L                   # consensus filter
)

obj <- loopQTL_normalize(
  obj,
  method  = "multiHiCcompare",      # or "BNBC"
  normal  = TRUE,                   # apply RINT after normalization
  BPPARAM = BiocParallel::MulticoreParam(workers = 8)
)

# Covariates: by default `run_QTL()` calls loopQTL_choose_K() (PCAForQTL elbow)
# and uses the chosen K phenotype PCs as MatrixEQTL's `cvrt`. To override,
# pass your own covariates matrix:
#   pcs <- loopQTL_choose_K(obj, K_method = "BE")$PCs
#   res <- run_QTL(obj, window = 0L, covariates = pcs)

results <- run_QTL(obj, window = 0L)
head(results)
#   loop_id   snp_id   chrom  pos    ref alt anchor   beta   se     t      p     FDR_cis  n_used  maf
```

## Run in parallel

The unit of parallel work is one `(sample, chromosome)` extraction in `loopQTL_setup()`. Dispatched via `BiocParallel::bplapply` — the backend is your choice via `BPPARAM`:

```r
# workstation:
BPPARAM = BiocParallel::MulticoreParam(workers = 16)

# SLURM cluster:
BPPARAM = BiocParallel::BatchtoolsParam(workers = 200, cluster = "slurm",
                                        template = "slurm.tmpl")
```
