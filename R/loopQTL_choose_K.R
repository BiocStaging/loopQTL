#' Choose the number of phenotype PCs to use as covariates for `run_QTL()`
#'
#' Wraps [PCAForQTL] (Zhou et al.) to pick how many phenotype PCs to include
#' as `cvrt` in MatrixEQTL. Default `K_method = "elbow"` (fast); `"BE"` is the
#' Buja-Eyuboglu permutation-based method (slower, more rigorous).
#'
#' `normal` controls whether RINT is applied to the phenotype before PCA.
#' Default `TRUE` matches `run_QTL()`'s default, so the PCs reflect the same
#' phenotype form that the regression will use. Pass `normal = FALSE` to PCA
#' the raw (BNBC/fastlo) phenotype instead.
#'
#' @param obj A [loopQTL-class] with `phenotype` populated by `loopQTL_normalize()`.
#' @param K_method `"elbow"` (default) or `"BE"`.
#' @param normal Apply RINT per loop before PCA. Default `TRUE`.
#' @param K_max Hard cap on the number of PCs returned. Default
#'   `length(obj@samples) - 2L`, which is the largest K that still leaves
#'   degrees of freedom for MatrixEQTL's intercept + genotype regressor.
#'   When the auto-selected K exceeds `K_max`, K is capped and a message
#'   is emitted. Pass a smaller value (e.g. `floor(n_samples / 4)`) to
#'   force a more conservative covariate set.
#' @param ... Forwarded to `PCAForQTL::runElbow` or `PCAForQTL::runBE`.
#'
#' @return A list:
#'   \describe{
#'     \item{`K`}{Chosen number of PCs (integer).}
#'     \item{`PCs`}{`samples x K` matrix of principal components (rownames = `samples`).}
#'     \item{`method`}{Which K-selection method was used.}
#'   }
#'
#' @export
loopQTL_choose_K <- function(obj,
                             K_method = c("elbow", "BE"),
                             normal   = TRUE,
                             K_max    = NULL,
                             ...) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    length(obj@phenotype) > 0L
  )
  K_method <- match.arg(K_method)

  if (is.null(K_max)) K_max <- max(1L, length(obj@samples) - 2L)
  K_max <- as.integer(K_max)

  Y <- obj@phenotype
  if (normal) Y <- .rint_matrix(Y)

  .pca_for_qtl(Y, samples = obj@samples, K_method = K_method,
               K_max = K_max, ...)
}

# Internal: PCA + K selection on a (loops x samples) matrix. Shared by
# loopQTL_choose_K() and run_QTL() (so we don't RINT twice when run_QTL
# already RINT'd the matrix it's passing).
.pca_for_qtl <- function(Y,
                         samples,
                         K_method = c("elbow", "BE"),
                         K_max    = NULL,
                         ...) {
  K_method <- match.arg(K_method)
  if (is.null(K_max)) K_max <- max(1L, length(samples) - 2L)
  K_max <- as.integer(K_max)

  keep_rows <- stats::complete.cases(Y)
  if (sum(keep_rows) < 2L)
    stop("Not enough non-NA loops to run PCA")
  Yt <- t(Y[keep_rows, , drop = FALSE])
  pca_res <- stats::prcomp(Yt, center = TRUE, scale. = FALSE)

  K_raw <- switch(
    K_method,
    elbow = PCAForQTL::runElbow(prcompResult = pca_res, ...),
    BE    = PCAForQTL::runBE(prcompResult = pca_res, ...)
  )
  K_raw <- as.integer(K_raw)
  if (is.na(K_raw) || K_raw < 1L) K_raw <- 1L
  K_raw <- min(K_raw, ncol(pca_res$x))

  K <- min(K_raw, K_max)
  if (K < K_raw) {
    message(sprintf(
      "PCAForQTL %s chose K = %d, capped to K_max = %d (n_samples - 2 leaves df for intercept + genotype).",
      K_method, K_raw, K))
  }

  PCs <- pca_res$x[, seq_len(K), drop = FALSE]
  rownames(PCs) <- samples
  colnames(PCs) <- paste0("PC", seq_len(K))

  list(K = K, PCs = PCs, method = K_method, K_raw = K_raw, K_max = K_max)
}
