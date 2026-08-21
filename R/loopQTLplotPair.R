#' Boxplot for a single SNP-loop pair, with QTL stats in the title
#'
#' Stratifies the loop's phenotype by genotype dosage (0/1/2, labeled with the
#' REF/ALT alleles) and adds a title showing snp_id, loop_id, REF/ALT, p-value,
#' and beta - drawn from the pair's row in `obj@results`.
#'
#' @param obj A [loopQTL-class] with `obj@results` populated by `runQTL()`.
#' @param loop_id Single loop_id to plot.
#' @param snp_id Single snp_id to plot.
#' @param normal Use RINT'd phenotype values (default `TRUE`, matching what
#'   the QTL regression actually used). Pass `FALSE` to plot pre-RINT values.
#' @param jitter Add a jittered point overlay (default `TRUE`).
#'
#' @return A `ggplot` object.
#'
#' @export
loopQTLplotPair <- function(obj, loop_id, snp_id,
                            normal = TRUE, jitter = TRUE) {
    stopifnot(
        methods::is(obj, "loopQTL"),
        is.character(loop_id), length(loop_id) == 1L,
        is.character(snp_id),  length(snp_id) == 1L
    )
    if (nrow(obj@results) == 0L) {
        .stopNoCall("obj@results is empty. Run runQTL(obj) first.")
    }

    ## Look up the (loop, SNP) row in obj@results for stats
    hit <- obj@results[obj@results$loop_id == loop_id &
        obj@results$snp_id == snp_id, , drop = FALSE]
    if (nrow(hit) == 0L) {
        .stopNoCall("(loop_id, snp_id) pair not found in obj@results")
    }
    hit <- hit[1L, , drop = FALSE]
    ref <- hit$ref
    alt <- hit$alt
    geno_labels <- c(
        paste0(ref, "/", ref),
        paste0(ref, "/", alt),
        paste0(alt, "/", alt)
    )

    ## Pull per-sample pair data
    d <- loopQTLpairData(obj, loop_id, snp_id, normal = normal)
    d$geno_factor <- factor(d$genotype, levels = c(0, 1, 2), labels = geno_labels)
    d <- d[!is.na(d$geno_factor) & !is.na(d$phenotype), , drop = FALSE]
    if (nrow(d) == 0L) {
        .stopNoCall("No samples have non-NA genotype and phenotype for this pair.")
    }

    ## Per-genotype sample counts for x-tick subtitles
    n_per <- as.integer(table(d$geno_factor))
    x_labels <- sprintf("%s\n(n=%d)", geno_labels, n_per)

    title_main <- sprintf("%s   |   %s", snp_id, loop_id)
    sub_main <- sprintf(
        "REF=%s   ALT=%s   |   p = %.3g   |   beta = %.3g",
        ref, alt, hit$p, hit$beta
    )
    y_lab <- if (normal) "Loop phenotype (RINT)" else "Loop phenotype (normalized)"

    p <- ggplot2::ggplot(d, ggplot2::aes(.data$geno_factor, .data$phenotype)) +
        ggplot2::geom_boxplot(
            fill = "lightgray", outlier.size = 0.8,
            width = 0.5
        ) +
        ggplot2::scale_x_discrete(labels = x_labels) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(
            title = title_main,
            subtitle = sub_main,
            x = "Genotype",
            y = y_lab
        ) +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 11),
            plot.subtitle = ggplot2::element_text(size = 10)
        )
    if (jitter) {
        p <- p + ggplot2::geom_jitter(
            width = 0.12, alpha = 0.6,
            size = 1.0, color = "steelblue"
        )
    }
    p
}
