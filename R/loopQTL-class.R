#' loopQTL S4 class
#'
#' Container for a loopQTL analysis. All matrices share the canonical sample
#' order in `samples`; this is an invariant every pipeline function preserves.
#'
#' @slot samples Canonical sample order.
#' @slot loops Consensus intra-chromosomal loops as a [InteractionSet::GInteractions].
#' @slot contacts Raw contact matrix (loops x samples), NA -> 0.
#' @slot cis Per-sample cis-contact totals, names match `samples`.
#' @slot vcf_path Path to the tabix-indexed VCF.
#' @slot resolution Analysis bin size (bp).
#' @slot phenotype Model-ready phenotype (loops x samples) after normalize();
#'   length-0 matrix before normalization.
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
    meta       = "list"
  ),
  prototype(
    samples    = character(0),
    contacts   = matrix(numeric(0), 0, 0),
    cis        = numeric(0),
    vcf_path   = NA_character_,
    resolution = NA_integer_,
    phenotype  = matrix(numeric(0), 0, 0),
    meta       = list()
  )
)

setValidity("loopQTL", function(object) {
  errs <- character(0)

  if (length(object@samples) == 0L)
    errs <- c(errs, "samples must be non-empty")

  if (anyDuplicated(object@samples))
    errs <- c(errs, "samples must be unique")

  if (!identical(colnames(object@contacts), object@samples))
    errs <- c(errs, "colnames(contacts) must equal samples")

  if (length(object@loops) != nrow(object@contacts))
    errs <- c(errs, "nrow(contacts) must equal length(loops)")

  if (length(object@cis) != length(object@samples) ||
      !identical(names(object@cis), object@samples))
    errs <- c(errs, "names(cis) must equal samples")

  if (length(object@phenotype) > 0L) {
    if (!identical(colnames(object@phenotype), object@samples))
      errs <- c(errs, "colnames(phenotype) must equal samples")
    if (nrow(object@phenotype) != nrow(object@contacts))
      errs <- c(errs, "phenotype and contacts must have the same number of rows")
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
  cat("  contacts  :", nrow(object@contacts), "x", ncol(object@contacts),
      "  range [", signif(min(object@contacts), 3), ",",
      signif(max(object@contacts), 3), "]\n")
  cat("  phenotype :",
      if (length(object@phenotype) > 0L)
        paste(nrow(object@phenotype), "x", ncol(object@phenotype))
      else "not yet computed", "\n")
  if (length(object@meta) > 0L) {
    cat("  meta      :\n")
    for (nm in names(object@meta))
      cat("    ", nm, "=", paste(object@meta[[nm]], collapse = ", "), "\n")
  }
  invisible(NULL)
})
