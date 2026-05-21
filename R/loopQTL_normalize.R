#' Normalize loopQTL contact counts to a model-ready phenotype
#'
#' Applies the originating-analysis transform in fixed order:
#' DESeq2 size factors -> log2(x+1) -> residualize on `log10(cis+1)` (and any
#' user covariates) -> RINT per loop. RINT is applied exactly once. Covariates
#' are "spent" here; [run_QTL()] must not re-include them.
#'
#' @param obj A [loopQTL-class] from [loopQTL_setup()].
#' @param covariates Optional matrix or data.frame of extra covariates
#'   (samples x covariates) to regress out alongside `log10(cis+1)`. Row names
#'   or row order must match `samples(obj)`.
#'
#' @return The input object with the `phenotype` slot filled. Zero-variance
#'   loops (residual sd == 0) are carried as all-NA rows; the count is recorded
#'   in `meta$n_zero_var`.
#'
#' @export
loopQTL_normalize <- function(obj, covariates = NULL) {
  stopifnot(methods::is(obj, "loopQTL"))
  samples_vec <- obj@samples
  contacts <- obj@contacts

  # ---- DESeq2 size factors + log2 ----
  sf <- DESeq2::estimateSizeFactorsForMatrix(contacts)
  norm_mat <- t(t(contacts) / sf)
  Y <- log2(norm_mat + 1)  # loops x samples

  # ---- Covariate design ----
  log_cis <- log10(obj@cis + 1)
  cov_mat <- data.frame(log_cis = log_cis, row.names = samples_vec)

  if (!is.null(covariates)) {
    extra <- as.data.frame(covariates, stringsAsFactors = FALSE)
    if (!is.null(rownames(extra))) {
      m <- match(samples_vec, rownames(extra))
      if (anyNA(m))
        stop("`covariates` is missing rows for samples: ",
             paste(samples_vec[is.na(m)], collapse = ", "))
      extra <- extra[m, , drop = FALSE]
    } else {
      if (nrow(extra) != length(samples_vec))
        stop("`covariates` has ", nrow(extra), " rows but there are ",
             length(samples_vec), " samples; provide rownames or matched order")
    }
    # drop constant / all-NA columns
    bad <- vapply(extra, function(x) all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0,
                  logical(1))
    if (any(bad)) {
      warning("Dropping non-informative covariate(s): ",
              paste(names(extra)[bad], collapse = ", "))
      extra <- extra[, !bad, drop = FALSE]
    }
    if (ncol(extra) > 0L) {
      rownames(extra) <- samples_vec
      cov_mat <- cbind(cov_mat, extra)
    }
  }

  X <- stats::model.matrix(~ ., data = cov_mat)

  if (qr(X)$rank < ncol(X))
    stop("Covariate design is rank-deficient; check for collinear covariates")

  # ---- QR residualization, all loops in one shot ----
  qrX <- qr(X)
  # fitted = X %*% qr.coef(qrX, t(Y))  -> p x n_loops, then transpose to align
  coefs <- qr.coef(qrX, t(Y))            # (p x n_loops)
  fitted <- t(X %*% coefs)               # (n_loops x n_samples)
  Y_resid <- Y - fitted

  # ---- Zero-variance detection ----
  row_sd <- apply(Y_resid, 1L, stats::sd)
  zero_var <- !is.finite(row_sd) | row_sd == 0
  n_zero_var <- sum(zero_var)

  # ---- RINT per loop ----
  phenotype <- matrix(NA_real_, nrow = nrow(Y_resid), ncol = ncol(Y_resid),
                      dimnames = list(NULL, samples_vec))
  if (any(!zero_var)) {
    good <- which(!zero_var)
    phenotype[good, ] <- t(apply(Y_resid[good, , drop = FALSE], 1L, .rint))
  }

  obj@phenotype <- phenotype
  obj@meta$n_zero_var <- n_zero_var
  obj@meta$normalize_order <- c("DESeq2", "log2", "residualize", "RINT")
  obj@meta$covariates <- colnames(cov_mat)
  methods::validObject(obj)
  obj
}
