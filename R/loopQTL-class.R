#' loopQTL S4 class
#'
#' Container for a loopQTL analysis. All matrices share the canonical sample
#' order in `samples`; this is an invariant every pipeline function preserves.
#'
#' @slot samples Canonical sample order.
#' @slot loops Consensus intra-chromosomal loops as a [InteractionSet::GInteractions].
#' @slot contacts Raw contact matrix (loops x samples), NA -> 0. Empty by
#'   default; populated only if you explicitly do a consensus extraction for QC.
#' @slot cis Per-sample cis-contact totals, names match `samples`. Empty
#'   unless the user supplied a `cis` column in the samplesheet.
#' @slot vcf_path Path to the tabix-indexed VCF.
#' @slot resolution Analysis bin size (bp).
#' @slot phenotype Model-ready phenotype (loops x samples) after normalize();
#'   length-0 matrix before normalization.
#' @slot results QTL result table written by `run_QTL()`. Empty data.frame
#'   before `run_QTL()` runs. Carries one row per tested (loop, SNP) pair
#'   with columns `loop_id`, `snp_id`, `chrom`, `pos`, `ref`, `alt`,
#'   `anchor`, `beta`, `se`, `t`, `p`, `FDR_cis`, `n_used`, `maf`. After
#'   `loopQTL_annotate()`, also carries `gene`, `peak` (if peak supplied),
#'   and `Type`.
#' @slot meta Free-form run metadata (counts of dropped loops, zero-var rows, ...).
#'
#' @name loopQTL-class
#' @rdname loopQTL-class
#' @exportClass loopQTL
setClass(
  "loopQTL",
  representation(
    samples    = "character",
    loops      = "GInteractions",
    contacts   = "matrix",
    cis        = "numeric",
    vcf_path   = "character",
    resolution = "integer",
    phenotype  = "matrix",
    results    = "data.frame",
    meta       = "list"
  ),
  prototype(
    samples    = character(0),
    contacts   = matrix(numeric(0), 0, 0),
    cis        = numeric(0),
    vcf_path   = NA_character_,
    resolution = NA_integer_,
    phenotype  = matrix(numeric(0), 0, 0),
    results    = data.frame(),
    meta       = list()
  )
)

setValidity("loopQTL", function(object) {
  errs <- character(0)

  if (length(object@samples) == 0L)
    errs <- c(errs, "samples must be non-empty")

  if (anyDuplicated(object@samples))
    errs <- c(errs, "samples must be unique")

  # contacts is optional (empty by default; populated only if the user
  # explicitly asks for raw extraction). When non-empty, columns must match
  # samples and rows must match loops.
  if (length(object@contacts) > 0L) {
    if (!identical(colnames(object@contacts), object@samples))
      errs <- c(errs, "colnames(contacts) must equal samples")
    if (length(object@loops) != nrow(object@contacts))
      errs <- c(errs, "nrow(contacts) must equal length(loops)")
  }

  # cis is optional. When present (length > 0), must align to samples.
  if (length(object@cis) > 0L) {
    if (length(object@cis) != length(object@samples) ||
        !identical(names(object@cis), object@samples))
      errs <- c(errs, "names(cis) must equal samples")
  }

  if (length(object@phenotype) > 0L) {
    if (!identical(colnames(object@phenotype), object@samples))
      errs <- c(errs, "colnames(phenotype) must equal samples")
    if (nrow(object@phenotype) != length(object@loops))
      errs <- c(errs, "phenotype must have one row per loop")
  }

  if (length(errs) == 0L) TRUE else errs
})

#' @rdname loopQTL-class
#' @param object a loopQTL object
#' @export
setMethod("show", "loopQTL", function(object) {
  cat("loopQTL object\n")
  cat("  samples   :", length(object@samples),
      "(", paste(utils::head(object@samples, 4), collapse = ", "),
      if (length(object@samples) > 4) ", ..." else "", ")\n")
  cat("  loops     :", length(object@loops), "consensus intra-chrom loops\n")
  cat("  resolution:", object@resolution, "bp\n")
  if (length(object@contacts) > 0L) {
    cat("  contacts  :", nrow(object@contacts), "x", ncol(object@contacts),
        "  range [", signif(min(object@contacts), 3), ",",
        signif(max(object@contacts), 3), "]\n")
  } else {
    cat("  contacts  : not extracted\n")
  }
  cat("  phenotype :",
      if (length(object@phenotype) > 0L)
        paste(nrow(object@phenotype), "x", ncol(object@phenotype))
      else "not yet computed", "\n")
  cat("  results   :",
      if (nrow(object@results) > 0L)
        paste(nrow(object@results), "(loop, SNP) tests")
      else "not yet computed", "\n")
  if (length(object@meta) > 0L) {
    cat("  meta      :\n")
    for (nm in names(object@meta))
      cat("    ", nm, "=", paste(object@meta[[nm]], collapse = ", "), "\n")
  }
  invisible(NULL)
})
