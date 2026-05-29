#' Choose the number of phenotype PCs to use as covariates for `run_QTL()`
#'
#' Wraps [PCAForQTL] (Zhou et al.) to pick how many phenotype PCs to include
#' as `cvrt` in MatrixEQTL. Default `K_method = "elbow"` (fast); `"BE"` is the
#' Buja-Eyuboglu permutation-based method (slower, more rigorous).
#'
#' @param obj A [loopQTL-class] with `phenotype` populated by `loopQTL_normalize()`.
#' @param K_method `"elbow"` (default) or `"BE"`.
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
                             ...) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    length(obj@phenotype) > 0L
  )
  K_method <- match.arg(K_method)

  # PCs of the samples are PCs of t(phenotype). Use only loops with finite values
  # (zero-var/NA rows from normalize were carried as NA — drop them here).
  Y <- obj@phenotype
  keep_rows <- stats::complete.cases(Y)
  if (sum(keep_rows) < 2L)
    stop("Not enough non-NA loops to run PCA")
  Yt <- t(Y[keep_rows, , drop = FALSE])
  pca_res <- stats::prcomp(Yt, center = TRUE, scale. = FALSE)

  K <- switch(
    K_method,
    elbow = PCAForQTL::runElbow(prcompResult = pca_res, ...),
    BE    = PCAForQTL::runBE(prcompResult = pca_res, ...)
  )
  K <- as.integer(K)
  if (is.na(K) || K < 1L) K <- 1L
  K <- min(K, ncol(pca_res$x))

  PCs <- pca_res$x[, seq_len(K), drop = FALSE]
  rownames(PCs) <- obj@samples
  colnames(PCs) <- paste0("PC", seq_len(K))

  list(K = K, PCs = PCs, method = K_method)
}
