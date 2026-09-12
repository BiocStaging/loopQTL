## ----echo = FALSE, message = FALSE--------------------------------------------
## `crop = NULL` resolves "The magick package is required to crop ..." error;
## see: https://stat.ethz.ch/pipermail/bioc-devel/2020-April/016656.html
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", crop = NULL)
options(tibble.print_min = 4, tibble.print_max = 4)

## ----eval = FALSE-------------------------------------------------------------
# ## Install BiocManager, which is required to install packages from Bioconductor
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# 
# BiocManager::install("loopQTL")

## ----eval = FALSE-------------------------------------------------------------
# ## Install prerequisite packages to install the development version from GitHub
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# if (!requireNamespace("remotes", quietly = TRUE)) {
#     install.packages("remotes")
# }
# 
# BiocManager::install(version = "devel")
# BiocManager::install("rhielab/loopQTL")

## ----setup, message = FALSE---------------------------------------------------
library(loopQTL)
library(ggplot2)

## -----------------------------------------------------------------------------
extdataDir <- system.file("extdata", package = "loopQTL")
vcfPath <- file.path(extdataDir, "example.vcf.gz")
samplesheetPath <- file.path(extdataDir, "example_samplesheet.tsv")
samples <- read.table(samplesheetPath, header = TRUE)

## Fix the input file paths to be absolute since we don't know what the working
## directory will be when the vignette is built. It is not necessary to do this
## during normal use of loopQTL if the working directory is set correctly.
samples$hic <- unlist(
    lapply(samples$hic, function(x) {
        file.path(extdataDir, x)
    })
)
samples$loops <- unlist(
    lapply(samples$loops, function(x) {
        file.path(extdataDir, x)
    })
)

## Display the format of the samplesheet
head(samples)

## -----------------------------------------------------------------------------
## Run loopQTLsetup
obj <- loopQTLsetup(
    samples    = samples,
    vcf        = vcfPath,
    resolution = 10000,
    minSamples = 2
)
obj

## -----------------------------------------------------------------------------
## Get the path to a temporary directory in which to save the output interact
## file
tempDirectory <- tempdir()

obj <- loopQTLnormalize(
    obj,
    BPPARAM     = BiocParallel::MulticoreParam(workers = 2),
    saveWindows = TRUE,
    windowsOut  = file.path(tempDirectory, "loop_snp_windows.rds"),
    window      = 2000000,
    snpWindow   = 10000
)

head(phenotype(obj)[, 1:5])

## -----------------------------------------------------------------------------
Kres <- loopQTLchooseK(obj, Kmethod = "elbow")
message("PCAForQTL chose K = ", Kres$K, " via ", Kres$method)
head(Kres$PCs)

## -----------------------------------------------------------------------------
obj <- runQTL(
    obj,
    window     = 10000,
    covariates = Kres$PCs
)
head(results(obj))

## -----------------------------------------------------------------------------
pair <- loopQTLpairData(
    obj,
    loopID = "L0000473_chr1_46050000_46130000",
    snpID  = "rs72677550"
)
head(pair)

## -----------------------------------------------------------------------------
peaksPath <- system.file("extdata", "example_annotation.bed", package = "loopQTL")
gtfPath <- system.file("extdata", "example.gtf.gz", package = "loopQTL")
obj <- loopQTLannotate(
    obj,
    peak            = peaksPath,
    peakAnnotation  = "Enhancer",
    TSSwindow       = 2000,
    removeTSS       = TRUE,
    gtf             = gtfPath
)
head(results(obj))

## -----------------------------------------------------------------------------
p <- loopQTLplotPair(
    obj,
    loopID = "L0000473_chr1_46050000_46130000",
    snpID  = "rs72677550"
)
p

