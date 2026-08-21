#' Per-pixel SNP-effect (beta) heatmap for one SNP
#'
#' For a SNP whose contact window was saved by
#' [loopQTLnormalize()] (with `save_windows = TRUE`), fits a per-cell
#' linear regression `contact_value ~ dosage [+ covariates]` and renders
#' the regression slope (`beta`) as a heatmap.
#'
#' Same units as the `beta` column in `obj@results` from [runQTL()] -
#' so with `normal = TRUE, covariates = "auto"` the cell at the loop's
#' (binA, binB) returns a beta that matches the QTL row directly (no
#' `r -> t` conversion needed). With `normal = FALSE, covariates = NULL`,
#' beta is the slope of the marginal regression on raw BNBC contact
#' values, matching the visualization-only style of Shi et al. 2025 but
#' on the slope scale instead of Pearson r.
#'
#' Optionally crop the view to a sub-window via `chrom`/`start`/`end` to
#' zoom on a specific region (e.g. just the loop arc).
#'
#' @param obj A [loopQTL-class] with `meta$windows_out` populated by
#'   `loopQTLnormalize(save_windows = TRUE)`.
#' @param snp_id Single snp_id present in the saved windows file.
#' @param chrom Optional chrom for verification. If supplied, must match
#'   the SNP's chrom in the windows file (helps catch a wrong snp_id).
#' @param start Optional start position in bp to crop the heatmap to.
#'   Snapped down to the nearest bin start.
#' @param end Optional end position in bp to crop the heatmap to. Snapped
#'   up to the nearest bin end.
#' @param windows Either `NULL` (default, loads from `obj@meta$windows_out`),
#'   a path to a windows `.rds`, or the loaded list. Pass the loaded list
#'   when calling the function in a loop over many SNPs to avoid reloading.
#' @param mirror Logical; mirror the upper-triangle beta values to the
#'   lower triangle so the heatmap is symmetric. Default `TRUE`.
#' @param minSamplesPerCell Minimum usable (non-NA dose + complete covariates) 
#'   samples required to fit the regression for a given cell. Cells with fewer
#'   are returned NA. Default 3.
#' @param normal Apply RINT (Blom-offset rank-inverse-normal transform) per
#'   cell across samples before fitting. Matches the `normal = TRUE`
#'   default in [runQTL()]. Default `FALSE` (regression on raw BNBC).
#' @param covariates Covariate handling for the regression:
#'   \describe{
#'     \item{`NULL`}{Default. Fit `cell ~ dosage` per cell (no covariates).}
#'     \item{`"auto"`}{Auto-pick K phenotype PCs via PCAForQTL on
#'       `obj@phenotype` (RINT'd if `normal = TRUE`), then fit
#'       `cell ~ dosage + PCs` per cell. Matches what [runQTL()] tested,
#'       so beta here equals the `beta` column in `obj@results`.}
#'     \item{matrix/data.frame}{User-supplied `samples x k` covariates;
#'       partial regression as above. Rownames or matched order required.}
#'   }
#' @param K_method `"elbow"` (default) or `"BE"`. Passed through to
#'   [loopQTLchooseK()] when `covariates = "auto"`. Ignored otherwise.
#' @param K_max Hard cap on the auto-chosen K when `covariates = "auto"`.
#'   Defaults to the number of samples minus 2 (same as [runQTL()]). Pass a
#'   smaller value (e.g. `floor(n_samples / 4)`) for a more conservative
#'   covariate set.
#' @param limit Symmetric color-scale limit on beta. Values beyond
#'   `+/-limit` are clipped (squished) to the deepest color. Default
#'   `NULL` = auto-set from `max(|beta|)` in the cropped view (beta is
#'   unbounded, unlike Pearson r, so there is no universally meaningful
#'   fixed limit). Pass an explicit numeric (e.g. `1.0`) to lock the
#'   scale across multiple SNPs for comparability.
#' @param palette RColorBrewer palette name for the divergent gradient.
#'   Default `"RdBu"` (reversed so + value = red, - value = blue).
#' @param metric What to render per cell:
#'   \describe{
#'     \item{`"beta_logp"`}{Default. Signed β weighted by significance:
#'       `beta * -log10(p)`. Highlights cells that are both large in
#'       effect and significant; cells with tiny p-values pop even if
#'       their raw β is modest.}
#'     \item{`"beta"`}{Raw regression slope β (same units as
#'       `obj@results$beta` when `normal = TRUE, covariates = "auto"`).}
#'   }
#' @param values Which cell values to run the regression on:
#'   \describe{
#'     \item{`"normalized"`}{Default. Post-BNBC values (`w$cells`). Uses the
#'       same normalization applied to `obj@phenotype`, so with
#'       `normal = TRUE, covariates = "auto"` the cell at the loop matches
#'       `obj@results$beta` exactly.}
#'     \item{`"cpm"`}{Counts per million: `w$cells_raw / (obj@cis / 1e6)`
#'       per sample. Depth-corrected but no cross-sample QN. Requires
#'       `obj@cis` to be populated (was set at setup).}
#'     \item{`"raw"`}{Raw straw counts (`w$cells_raw`) - no depth correction,
#'       no QN. β is on the raw-count scale and is NOT comparable to
#'       `obj@results$beta`.}
#'   }
#' @param axis_unit Unit for the graph axes - `"Mb"` (default) or `"bp"`.
#'
#' @return A `ggplot` object.
#'
#' @export
loopQTLplotSnpWindow <- function(obj, snp_id,
                                 chrom = NULL,
                                 start = NULL,
                                 end = NULL,
                                 windows = NULL,
                                 mirror = TRUE,
                                 minSamplesPerCell = 3L,
                                 normal = FALSE,
                                 covariates = NULL,
                                 K_method = c("elbow", "BE"),
                                 K_max = NULL,
                                 limit = NULL,
                                 palette = "RdBu",
                                 metric = c("beta_logp", "beta"),
                                 values = c("normalized", "cpm", "raw"),
                                 axis_unit = c("Mb", "bp")) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(snp_id), length(snp_id) == 1L, nzchar(snp_id),
        is.logical(mirror), length(mirror) == 1L,
        is.logical(normal), length(normal) == 1L, !is.na(normal),
        is.character(palette), length(palette) == 1L
    )
    if (!is.null(limit)) {
        stopifnot(is.numeric(limit), length(limit) == 1L, limit > 0)
    }
    K_method <- match.arg(K_method)
    metric <- match.arg(metric)
    values <- match.arg(values)
    axis_unit <- match.arg(axis_unit)

    w <- .resolve_windows(obj, windows)

    # 2. Look up the SNP row
    i <- match(snp_id, w$snp_index$snp_id)
    if (is.na(i)) {
        stop("snp_id not found in windows file: ", snp_id)
    }
    row <- w$snp_index[i, , drop = FALSE]
    cols <- row$col_indices[[1L]]
    if (length(cols) == 0L) {
        stop(
            "No cells were saved for snp_id = ", snp_id,
            " (likely on a chrom that had an empty window slice)."
        )
    }
    if (!is.null(chrom)) {
        if (length(chrom) != 1L || !is.character(chrom)) {
            stop("`chrom` must be a single character string")
        }
        if (!identical(chrom, row$chrom)) {
            stop(
                "Supplied chrom = '", chrom, "' does not match the SNP's chrom ('",
                row$chrom, "') in the windows file"
            )
        }
    }

    res <- as.integer(w$resolution)
    bin <- as.integer(row$binStart)
    binEnd <- as.integer(row$binEnd)

    # 3. Optionally crop to user's [start, end]
    if (!is.null(start) || !is.null(end)) {
        if (is.null(start)) start <- binStart
        if (is.null(end)) end <- binEnd + res
        if (!is.numeric(start) || !is.numeric(end) ||
            length(start) != 1L || length(end) != 1L || start >= end) {
            stop("`start` and `end` must be single numerics with start < end")
        }
        crop_lo_req <- as.integer(floor(start / res) * res)
        crop_hi_req <- as.integer(floor((end - 1L) / res) * res)
        crop_lo <- max(crop_lo_req, binStart)
        crop_hi <- min(crop_hi_req, binEnd)
        if (crop_lo > crop_hi) {
            stop(
                "Crop range [", start, ", ", end, "] does not overlap the saved ",
                "window [", binStart, ", ", binEnd + res, "] for SNP ", snp_id
            )
        }
        if (crop_lo != crop_lo_req || crop_hi != crop_hi_req) {
            warning(sprintf(
                paste0(
                    "Requested [%d, %d] was clipped to [%d, %d] (the saved window ",
                    "for %s covers [%d, %d]; widen `window` in loopQTLnormalize() ",
                    "to enlarge it)."
                ),
                as.integer(start), as.integer(end),
                crop_lo, crop_hi + res, snp_id, binStart, binEnd + res
            ))
        }
    } else {
        crop_lo <- binStart
        crop_hi <- binEnd
    }

    # 4. Pull genotype dosage for this SNP
    dose <- .read_snp_dosage(obj@vcf_path, row$chrom, row$pos,
        samples_vec = w$samples
    )

    # 5. Per-cell regression slope (beta)
    # Three modes, matching runQTL()'s semantics:
    #   normal = FALSE, covariates = NULL          -> marginal beta on raw BNBC
    #   normal = TRUE,  covariates = NULL          -> marginal beta on RINT'd cells
    #   normal = TRUE,  covariates = "auto"/matrix -> partial beta matching runQTL()
    # By Frisch-Waugh-Lovell, the partial slope equals the slope from
    # lm(cell ~ dose + cvrt), so we residualize both cell-value and dose
    # against [1, cvrt] then take slope = sum(Y_res * d_res) / sum(d_res^2).
    # Select the cell values matrix based on the `values` arg
    
    if (values %in% c("raw", "cpm") && is.null(w$cells_raw)) {
        .stopNoCall(
            "This loopQTL object does not contain raw counts. To run ",
            "loopQTLplotSnpWindow(), please re-run loopQTLnormalize() ",
            "with save_windows set to TRUE."
        )
    }
    
    cells_mat <- switch(values,
        normalized = w$cells,
        raw = w$cells_raw,
        cpm = {
            if (length(obj@cis) == 0L) {
                stop(
                    "values = 'cpm' requires obj@cis (per-sample cis depth). ",
                    "Set the `cis` column in the samplesheet at loopQTLsetup() or ",
                    "let setup auto-compute it from .hic."
                )
            }
            lib_M <- obj@cis[w$samples] / 1e6
            sweep(w$cells_raw, 1, lib_M, "/")
        }
    )
    cells_for_snp <- cells_mat[, cols, drop = FALSE]
    rownames(cells_for_snp) <- w$samples

    Y <- cells_for_snp
    if (normal) {
        Y <- apply(Y, 2L, .rint)
        if (!is.matrix(Y)) Y <- as.matrix(Y)
        if (ncol(Y) != length(cols)) Y <- t(Y)
    }

    cvrt_mode <- "none"
    cvrt <- NULL
    if (identical(covariates, "auto")) {
        cvrt_mode <- "auto"
        if (length(obj@phenotype) == 0L) {
            stop("covariates = 'auto' needs obj@phenotype; run loopQTLnormalize() first")
        }
        Y_for_pca <- if (isTRUE(normal)) .rint_matrix(obj@phenotype) else obj@phenotype
        K_res <- .pca_for_qtl(Y_for_pca,
            samples = obj@samples,
            K_method = K_method, K_max = K_max
        )
        cvrt <- K_res$PCs[w$samples, , drop = FALSE]
    } else if (!is.null(covariates)) {
        cvrt_mode <- "user"
        cv_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
        if (!is.null(rownames(cv_df))) {
            mm <- match(w$samples, rownames(cv_df))
            if (anyNA(mm)) {
                stop(
                    "`covariates` is missing rows for samples: ",
                    paste(w$samples[is.na(mm)], collapse = ", ")
                )
            }
            cv_df <- cv_df[mm, , drop = FALSE]
        } else if (nrow(cv_df) != length(w$samples)) {
            stop(
                "`covariates` has ", nrow(cv_df), " rows but ",
                length(w$samples), " samples; provide rownames or matched order"
            )
        }
        cvrt <- as.matrix(cv_df)
        rownames(cvrt) <- w$samples
    }

    usable <- !is.na(dose)
    if (!is.null(cvrt)) usable <- usable & stats::complete.cases(cvrt)
    n_usable <- sum(usable)

    beta_vec <- rep(NA_real_, ncol(Y))
    p_vec <- rep(NA_real_, ncol(Y))
    if (n_usable >= as.integer(minSamplesPerCell)) {
        Y_u <- Y[usable, , drop = FALSE]
        d_u <- dose[usable]
        if (is.null(cvrt)) {
            # Marginal: center then slope = sum(Yc * dc) / sum(dc^2)
            Y_r <- sweep(Y_u, 2L, colMeans(Y_u, na.rm = TRUE), FUN = "-")
            d_r <- d_u - mean(d_u)
            df <- n_usable - 2L
        } else {
            # Partial via FWL: residualize Y and d against [1, cvrt[usable, ]]
            X <- cbind(1, cvrt[usable, , drop = FALSE])
            qrX <- qr(X)
            Y_r <- Y_u - X %*% qr.coef(qrX, Y_u)
            d_r <- as.numeric(d_u - X %*% qr.coef(qrX, d_u))
            df <- n_usable - ncol(X) - 1L
        }
        denom_d <- sum(d_r^2)
        if (denom_d > 0 && df > 0) {
            beta_vec <- as.numeric(crossprod(Y_r, d_r) / denom_d)
            # RSS_j = ||Y_r_j||^2 - beta_j^2 * denom_d  (from partial regression)
            rss <- colSums(Y_r^2) - beta_vec^2 * denom_d
            rss <- pmax(rss, 0) # numerical safety
            sigma2 <- rss / df
            se_vec <- sqrt(sigma2 / denom_d)
            t_vec <- beta_vec / se_vec
            p_vec <- 2 * stats::pt(-abs(t_vec), df = df)
        }
        beta_vec[!is.finite(beta_vec)] <- NA_real_
        p_vec[!is.finite(p_vec)] <- NA_real_
    }

    # Build the display value per cell based on `metric`
    if (metric == "beta_logp") {
        # Guard against p==0 (perfect fit); cap the -log10 to a plottable ceiling
        p_safe <- pmax(p_vec, .Machine$double.xmin)
        display_v <- beta_vec * (-log10(p_safe))
        display_v[!is.finite(display_v)] <- NA_real_
    } else {
        display_v <- beta_vec
    }

    # 6. Lay r values onto the bins x bins grid
    coord <- w$cell_coords[cols, , drop = FALSE]
    bins <- seq.int(crop_lo, crop_hi, by = res)
    nb <- length(bins)
    in_crop <- coord$bin_i >= crop_lo & coord$bin_i <= crop_hi &
        coord$bin_j >= crop_lo & coord$bin_j <= crop_hi
    ri <- match(coord$bin_i[in_crop], bins)
    rj <- match(coord$bin_j[in_crop], bins)
    vv <- display_v[in_crop]

    # Shift tile centers by res/2 so each tile visually spans [bin_start,
    # bin_start + res] and axis ticks at bin_start align with the LEFT
    # edge of the tile -- matching IGV / cooler / Juicebox bp conventions.
    half <- res / 2
    long <- data.frame(
        bin_i = c(bins[ri] + half, if (mirror) bins[rj] + half),
        bin_j = c(bins[rj] + half, if (mirror) bins[ri] + half),
        value = c(vv, if (mirror) vv),
        stringsAsFactors = FALSE
    )
    long <- long[!duplicated(paste(long$bin_i, long$bin_j, sep = "_")), ,
        drop = FALSE
    ]

    # 7. Color scale (RdBu divergent on `value`)  # value scale depends on `metric`; both are unbounded. Default limit is
    # data-driven; pass a fixed `limit` for cross-SNP comparability.
    if (is.null(limit)) {
        vfin <- long$value[is.finite(long$value)]
        eff_limit <- if (length(vfin) > 0L) max(abs(vfin)) else 1
        if (!is.finite(eff_limit) || eff_limit <= 0) eff_limit <- 1
    } else {
        eff_limit <- limit
    }

    # 8. Plot
    axis_format <- if (axis_unit == "Mb") {
        function(x) sprintf("%.2f", x / 1e6)
    } else {
        function(x) format(x, big.mark = ",", scientific = FALSE)
    }
    axis_label_suffix <- if (axis_unit == "Mb") "Mb" else "bp"

    mode_str <- sprintf(
        "%s%s",
        if (cvrt_mode == "none") "marginal β" else "partial β",
        if (isTRUE(normal)) ", RINT" else ""
    )
    if (cvrt_mode == "auto") mode_str <- paste0(mode_str, sprintf(" | K=%d auto PCs", ncol(cvrt)))
    if (cvrt_mode == "user") mode_str <- paste0(mode_str, sprintf(" | %d user cvrt", ncol(cvrt)))
    mode_str <- paste0(mode_str, " | ", values, " values")

    metric_label <- switch(metric,
        beta_logp = expression(beta %*% -log[10](p)),
        beta      = expression(beta)
    )
    metric_title_str <- switch(metric,
        beta_logp = "β × -log10(p)",
        beta      = "β"
    )

    title <- sprintf(
        "Per-pixel %s (contact ~ dosage) | %s",
        metric_title_str, snp_id
    )
    subtitle <- sprintf(
        "%s:%d  REF=%s ALT=%s | window %s bp | %d samples (%d usable) | %s | limits ±%.3g",
        row$chrom, row$pos, row$ref, row$alt,
        format(crop_hi - crop_lo + res, big.mark = ","),
        length(w$samples), n_usable, mode_str, eff_limit
    )

    diag_lo <- crop_lo
    diag_hi <- crop_hi + res

    ggplot2::ggplot(
        long,
        ggplot2::aes(.data$bin_i, .data$bin_j, fill = .data$value)
    ) +
        ggplot2::geom_raster(interpolate = FALSE) +
        ggplot2::annotate("segment",
            x = diag_lo, y = diag_lo,
            xend = diag_hi, yend = diag_hi,
            color = "black", linewidth = 0.5, alpha = 0.8
        ) +
        ggplot2::scale_fill_distiller(
            palette   = palette,
            direction = -1, # so + beta = red, - beta = blue
            limits    = c(-eff_limit, eff_limit),
            oob       = scales::squish,
            na.value  = "grey95"
        ) +
        ggplot2::scale_x_continuous(labels = axis_format, expand = c(0, 0)) +
        ggplot2::scale_y_reverse(labels = axis_format, expand = c(0, 0)) +
        ggplot2::coord_fixed() +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = sprintf("bin_i (%s, %s)", axis_label_suffix, row$chrom),
            y = sprintf("bin_j (%s, %s)", axis_label_suffix, row$chrom),
            fill = metric_label
        ) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(
            plot.subtitle    = ggplot2::element_text(size = 9),
            panel.grid       = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            axis.line        = ggplot2::element_line(linewidth = 0.3)
        )
}

# Internal helpers

# Resolve the `windows` arg into a loaded list, with the necessary structure.
.resolve_windows <- function(obj, windows) {
    ## If the windows argument is null, look inside the loopQTL object
    if (is.null(windows)) {
        windows <- obj@meta$windows_out
        if (is.null(windows)) {
            .stopNoCall(
                "The input loopQTL object does not contain a windows file ",
                "path. Please run loopQTLnormalize() with save_windows set ",
                "to TRUE and windows_out set to an output file path, or ",
                "specify a path to a windows file via the windows argument."
            )
        }
    }
    if (is.character(windows) && length(windows) == 1L) {
        if (!file.exists(windows)) {
            stop("Windows file does not exist: ", windows)
        }
        return(readRDS(windows))
    }
    if (is.list(windows) &&
        all(c("cells", "cell_coords", "snp_index", "samples", "resolution") %in%
            names(windows))) {
        return(windows)
    }
    stop(
        "`windows` must be NULL, a path to a windows .rds, or a loaded list ",
        "with cells/cell_coords/snp_index/samples/resolution fields."
    )
}

# Targeted single-SNP dosage read from a tabix VCF. Returns a numeric vector
# of length(samples_vec) in that order; NA where missing.
.read_snp_dosage <- function(vcf_path, chrom, pos, samples_vec) {
    tab <- Rsamtools::TabixFile(vcf_path)
    vcf_chroms <- Rsamtools::headerTabix(tab)$seqnames
    vcf_norm <- .make_chrom_normalizer(vcf_chroms)
    gr <- GenomicRanges::GRanges(
        seqnames = vcf_norm(chrom),
        ranges   = IRanges::IRanges(start = as.integer(pos), end = as.integer(pos))
    )
    vcf <- VariantAnnotation::readVcf(
        vcf_path,
        genome = "unknown",
        param = VariantAnnotation::ScanVcfParam(
            which = gr, geno = "GT",
            info = NA_character_
        )
    )
    if (length(vcf) == 0L) {
        stop("VCF has no record at ", chrom, ":", pos)
    }
    # If multi-allelic, take the first row (matching how .gt_to_dosage handles it)
    gt <- VariantAnnotation::geno(vcf)$GT
    m <- match(samples_vec, colnames(gt))
    if (anyNA(m)) {
        stop(
            "VCF missing samples present in windows file: ",
            paste(samples_vec[is.na(m)], collapse = ", ")
        )
    }
    gt <- gt[1L, m, drop = FALSE]
    dose <- .gt_to_dosage(matrix(as.character(gt), nrow = 1L))[1L, ]
    setNames(as.numeric(dose), samples_vec)
}
