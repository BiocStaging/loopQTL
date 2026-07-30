#' Locus plot: RefGene + SNP + genotype-stratified loop arcs
#'
#' Stacks three tracks over a shared genomic-coordinate axis:
#'
#' 1. **Gene track** - genes from `gtf` overlapping `[start, end]`, drawn as
#'    strand-oriented rectangles with labels.
#' 2. **SNP track** - vertical marker at the SNP's position.
#' 3. **Loop arcs** - one panel per genotype class (0, 1, 2). For each
#'    significant loop anchored by the SNP (per `obj@results`), a Bezier
#'    arc connects the two anchor midpoints, colored by the loop's mean
#'    contact strength across samples in that genotype class.
#'
#' @param obj A [loopQTL-class] with `obj@results` populated by [runQTL()].
#' @param snp_id Single snp_id present in `obj@results`.
#' @param chrom Chromosome (must equal the SNP's chrom in `obj@results`).
#' @param start,end Genomic window in bp for the x-axis.
#' @param gtf Path to a GTF/GFF file (read via `rtracklayer::import()`).
#' @param p_threshold Nominal p-value cutoff for "significant" loops. Only
#'   loops with `obj@results$p <= p_threshold` (for this SNP) are drawn.
#'   Default `0.05`.
#' @param values Which contact-strength matrix to color arcs by:
#'   \describe{
#'     \item{`"normalized"`}{Default. Post-BNBC values from `obj@phenotype`.
#'       Matches what `runQTL()` tested.}
#'     \item{`"cpm"`}{`w$cells_raw / (obj@cis / 1e6)` per sample.
#'       Requires `save_windows` output; pass `windows` argument.}
#'     \item{`"raw"`}{Integer straw counts from `w$cells_raw`. Depth-
#'       confounded across samples; use with care. Requires `windows`.}
#'   }
#' @param windows Optional path/loaded list from `loopQTLnormalize(
#'   save_windows = TRUE)`. Required when `values %in% c("cpm", "raw")`.
#' @param strong_color Color for the strongest signal on the arc gradient.
#'   Default `"#FF7917"` (orange).
#' @param weak_color Color for the weakest signal. Default `"grey90"`.
#' @param limit Optional numeric cap on the color scale (values above are
#'   squished to `strong_color`). Default `NULL` = data-driven max.
#' @param gene_biotypes Character vector of `gene_biotype` values to keep
#'   from the GTF. Default `"protein_coding"`; pass `NULL` to keep all.
#' @param axis_unit `"Mb"` (default) or `"bp"` for x-axis tick labels.
#'
#' @return A `patchwork` composite `ggplot` (gene / SNP / arc panels
#'   stacked). Save with `ggplot2::ggsave()`.
#'
#' @export
loopQTLplotLocus <- function(obj, snp_id,
                                chrom, start, end,
                                gtf,
                                p_threshold  = 0.05,
                                values       = c("normalized", "cpm", "raw"),
                                windows      = NULL,
                                strong_color = "#FF7917",
                                weak_color   = "grey90",
                                limit        = NULL,
                                gene_biotypes = "protein_coding",
                                axis_unit    = c("Mb", "bp")) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    is.character(snp_id), length(snp_id) == 1L, nzchar(snp_id),
    is.character(chrom),  length(chrom)  == 1L,
    is.numeric(start), length(start) == 1L,
    is.numeric(end),   length(end)   == 1L, start < end,
    is.character(gtf), length(gtf) == 1L, file.exists(gtf),
    is.numeric(p_threshold), length(p_threshold) == 1L,
    p_threshold > 0, p_threshold <= 1
  )
  values    <- match.arg(values)
  axis_unit <- match.arg(axis_unit)
  for (pkg in c("ggplot2", "patchwork", "rtracklayer", "scales"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " is required for loopQTLplotLocus(); install it first.")
  if (nrow(obj@results) == 0L)
    stop("obj@results is empty. Run runQTL(obj) first.")

  start <- as.integer(start); end <- as.integer(end)

  # 1. Filter obj@results to (this SNP, p <= threshold) and to loops in the region
  resultsForThisSNP <- obj@results[obj@results$snp_id == snp_id, , drop = FALSE]
  if (nrow(resultsForThisSNP) == 0L)
    stop("SNP ", snp_id, " not found in obj@results")
  if (!identical(unique(as.character(resultsForThisSNP$chrom)), chrom))
    stop("Supplied chrom = '", chrom, "' does not match the SNP's chrom ('",
         resultsForThisSNP$chrom[1], "') in obj@results")
  snp_pos <- resultsForThisSNP$pos[1]
  snp_ref <- resultsForThisSNP$ref[1]
  snp_alt <- resultsForThisSNP$alt[1]

  significantResults <- resultsForThisSNP[resultsForThisSNP$p <= p_threshold & is.finite(resultsForThisSNP$p), ,
                     drop = FALSE]
  if (nrow(significantResults) == 0L)
    stop("No significant (p <= ", p_threshold, ") loops for SNP ", snp_id)

  # Look up anchor coords for each significant loop; filter to region
  loop_ids_all <- S4Vectors::mcols(obj@loops)$loop_id
  ancA_all <- InteractionSet::anchors(obj@loops, "first")
  ancB_all <- InteractionSet::anchors(obj@loops, "second")
  li <- match(significantResults$loop_id, loop_ids_all)
  loops_df <- data.frame(
    loop_id = significantResults$loop_id,
    chrom   = as.character(GenomicRanges::seqnames(ancA_all[li])),
    ancA_start = GenomicRanges::start(ancA_all[li]),
    ancA_end   = GenomicRanges::end(ancA_all[li]),
    ancB_start = GenomicRanges::start(ancB_all[li]),
    ancB_end   = GenomicRanges::end(ancB_all[li]),
    beta = significantResults$beta,
    p    = significantResults$p,
    row_i = li,
    stringsAsFactors = FALSE
  )
  in_region <- loops_df$chrom == chrom &
               pmin(loops_df$ancA_start, loops_df$ancB_start) >= start &
               pmax(loops_df$ancA_end,   loops_df$ancB_end)   <= end
  loops_df <- loops_df[in_region, , drop = FALSE]
  if (nrow(loops_df) == 0L)
    stop("No significant loops fully inside [", start, ", ", end,
         "] for SNP ", snp_id)
  loops_df$ancA_mid <- (loops_df$ancA_start + loops_df$ancA_end) / 2
  loops_df$ancB_mid <- (loops_df$ancB_start + loops_df$ancB_end) / 2

  # 2. Per-sample dosage for this SNP
  dose <- .read_snp_dosage(obj@vcf_path, chrom, snp_pos,
                           samples_vec = obj@samples)
  # Split samples into genotype classes 0/1/2 (drop NA)
  geno_groups <- split(names(dose)[!is.na(dose)],
                       factor(dose[!is.na(dose)], levels = c(0, 1, 2)))

  # 3. For each loop, get per-sample contact values, then mean per genotype
  values_mat <- .locus_values_matrix(obj, loops_df, values, windows)

  mean_per_geno <- vapply(names(geno_groups), function(g) {
    ss <- geno_groups[[g]]
    if (length(ss) == 0L) return(rep(NA_real_, nrow(loops_df)))
    rowMeans(values_mat[, ss, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(loops_df)))
  colnames(mean_per_geno) <- names(geno_groups)

  n_per <- vapply(geno_groups, length, integer(1L))
  geno_labels <- c(paste0(snp_ref, "/", snp_ref),
                   paste0(snp_ref, "/", snp_alt),
                   paste0(snp_alt, "/", snp_alt))
  panel_titles <- sprintf("%s (n=%d)", geno_labels, n_per)

  # 4. Color scale: use ACTUAL data range so small differences show up.
  # BNBC values can be small (or negative); a fixed 0..max scale often
  # squashes real differences into the pale end of the gradient.
  vfin <- mean_per_geno[is.finite(mean_per_geno)]
  if (!is.null(limit)) {
    eff_lo <- 0; eff_hi <- limit
  } else if (length(vfin) > 0L) {
    eff_lo <- min(vfin); eff_hi <- max(vfin)
    if (eff_hi - eff_lo < .Machine$double.eps) {
      eff_lo <- eff_lo - 0.5; eff_hi <- eff_hi + 0.5
    }
  } else {
    eff_lo <- 0; eff_hi <- 1
  }

  # 5. Read GTF genes + exons in region (IGV-style: exons as blocks,
  #    introns as thin lines, strand arrows)
  gene_df <- .load_genes_in_region(gtf, chrom, start, end, gene_biotypes)
  exon_df <- .load_exons_in_region(gtf, chrom, start, end, gene_df$gene)

  # 6. Axis formatter
  ax_fmt <- if (axis_unit == "Mb")
    function(x) sprintf("%.2f", x / 1e6) else
    function(x) format(x, big.mark = ",", scientific = FALSE)

  # 7. Build the four ggplot panels
  p_gene <- .build_gene_track(gene_df, exon_df, start, end, chrom, ax_fmt)
  p_snp  <- .build_snp_track(snp_pos, snp_id, start, end, ax_fmt)

  arc_panels <- lapply(seq_along(geno_groups), function(g_i) {
    g_lab <- names(geno_groups)[g_i]
    arc_df <- data.frame(
      loop_id = loops_df$loop_id,
      x1      = loops_df$ancA_mid,
      x2      = loops_df$ancB_mid,
      value   = mean_per_geno[, g_lab]
    )
    .build_arc_panel(arc_df, start, end, ax_fmt,
                     eff_lo, eff_hi, strong_color, weak_color,
                     panel_titles[g_i],
                     show_x = (g_i == length(geno_groups)))
  })

  # 8. Stack with patchwork: gene | SNP | genotype 0 | 1 | 2
  #    Heights: gene 1.2, SNP 0.4, arcs each 1.2
  #    guides = "collect" consolidates all three arc legends into one.
  heights <- c(1.2, 0.4, rep(1.2, length(arc_panels)))
  wrap    <- patchwork::wrap_plots(
    c(list(p_gene, p_snp), arc_panels),
    ncol = 1,
    heights = heights
  ) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = sprintf("%s | %s:%d-%d",
                         snp_id, chrom, start, end),
      subtitle = sprintf("REF=%s  ALT=%s  |  %d significant loop(s) (p <= %.3g)  |  %s values",
                         snp_ref, snp_alt, nrow(loops_df),
                         p_threshold, values)
    )
  wrap
}


# Internal helpers

# Pull per-sample contact values for the given loops based on `values`.
# Returns a (n_loops x n_samples) matrix in obj@samples order.
.locus_values_matrix <- function(obj, loops_df, values, windows) {
  loop_ids_all <- S4Vectors::mcols(obj@loops)$loop_id
  li <- match(loops_df$loop_id, loop_ids_all)
  if (values == "normalized") {
    if (length(obj@phenotype) == 0L)
      stop("values = 'normalized' needs obj@phenotype; run loopQTLnormalize() first.")
    return(obj@phenotype[li, obj@samples, drop = FALSE])
  }
  # cpm / raw modes need windows data
  w <- .resolve_windows(obj, windows)
  if (is.null(w$cells_raw))
    stop("values = '", values, "' needs `cells_raw` in the windows file; ",
         "re-run loopQTLnormalize(save_windows = TRUE).")
  res <- obj@resolution
  n_loops <- nrow(loops_df)
  out <- matrix(0, nrow = n_loops, ncol = length(obj@samples),
                dimnames = list(loops_df$loop_id, obj@samples))
  for (i in seq_len(n_loops)) {
    binA <- .bin_anchor(loops_df$ancA_start[i], loops_df$ancA_end[i], res)
    binB <- .bin_anchor(loops_df$ancB_start[i], loops_df$ancB_end[i], res)
    lo <- min(binA, binB); hi <- max(binA, binB)
    hit <- which(w$cell_coords$chrom == loops_df$chrom[i] &
                 w$cell_coords$bin_i == lo &
                 w$cell_coords$bin_j == hi)
    if (length(hit) == 0L) next
    raw_vals <- w$cells_raw[obj@samples, hit[1]]
    if (values == "cpm") {
      if (length(obj@cis) == 0L)
        stop("values = 'cpm' needs obj@cis (per-sample cis depth).")
      raw_vals <- raw_vals / (obj@cis[obj@samples] / 1e6)
    }
    out[i, ] <- raw_vals
  }
  out
}

# Read genes overlapping [start, end] on chrom from a GTF. Handles chr vs no-chr.
.load_genes_in_region <- function(gtf, chrom, start, end, biotypes) {
  gr <- rtracklayer::import(gtf)
  type_col <- if ("type" %in% names(S4Vectors::mcols(gr))) "type"
              else if ("feature" %in% names(S4Vectors::mcols(gr))) "feature"
              else stop("GTF has neither `type` nor `feature` column")
  genes <- gr[S4Vectors::mcols(gr)[[type_col]] == "gene"]
  if (length(genes) == 0L) stop("No `gene` entries in GTF")

  # Chrom naming
  gtf_chroms <- unique(as.character(GenomicRanges::seqnames(genes)))
  chrom_q <- if (any(grepl("^chr", gtf_chroms)))
    (if (grepl("^chr", chrom)) chrom else paste0("chr", chrom))
  else sub("^chr", "", chrom)

  keep <- as.character(GenomicRanges::seqnames(genes)) == chrom_q &
          GenomicRanges::end(genes)   >= start &
          GenomicRanges::start(genes) <= end
  genes <- genes[keep]
  if (!is.null(biotypes) && "gene_biotype" %in% names(S4Vectors::mcols(genes))) {
    genes <- genes[as.character(S4Vectors::mcols(genes)$gene_biotype) %in% biotypes]
  }

  name_col <- if ("gene_name" %in% names(S4Vectors::mcols(genes))) "gene_name"
              else if ("gene_id" %in% names(S4Vectors::mcols(genes))) "gene_id"
              else NULL
  data.frame(
    gene   = if (is.null(name_col)) NA_character_
             else as.character(S4Vectors::mcols(genes)[[name_col]]),
    start  = pmax(start, GenomicRanges::start(genes)),
    end    = pmin(end,   GenomicRanges::end(genes)),
    strand = as.character(GenomicRanges::strand(genes)),
    stringsAsFactors = FALSE
  )
}

# Read exons of the given gene names in the region.
.load_exons_in_region <- function(gtf, chrom, start, end, gene_names) {
  if (length(gene_names) == 0L)
    return(data.frame(gene = character(0), start = integer(0),
                      end = integer(0), strand = character(0),
                      stringsAsFactors = FALSE))
  gr <- rtracklayer::import(gtf)
  type_col <- if ("type" %in% names(S4Vectors::mcols(gr))) "type"
              else if ("feature" %in% names(S4Vectors::mcols(gr))) "feature"
              else return(data.frame(gene = character(0), start = integer(0),
                                     end = integer(0), strand = character(0),
                                     stringsAsFactors = FALSE))
  exons <- gr[S4Vectors::mcols(gr)[[type_col]] == "exon"]
  if (length(exons) == 0L)
    return(data.frame(gene = character(0), start = integer(0),
                      end = integer(0), strand = character(0),
                      stringsAsFactors = FALSE))
  gtf_chroms <- unique(as.character(GenomicRanges::seqnames(exons)))
  chrom_q <- if (any(grepl("^chr", gtf_chroms)))
    (if (grepl("^chr", chrom)) chrom else paste0("chr", chrom))
  else sub("^chr", "", chrom)
  name_col <- if ("gene_name" %in% names(S4Vectors::mcols(exons))) "gene_name"
              else if ("gene_id" %in% names(S4Vectors::mcols(exons))) "gene_id"
              else NULL
  if (is.null(name_col))
    return(data.frame(gene = character(0), start = integer(0),
                      end = integer(0), strand = character(0),
                      stringsAsFactors = FALSE))
  keep <- as.character(GenomicRanges::seqnames(exons)) == chrom_q &
          GenomicRanges::end(exons)   >= start &
          GenomicRanges::start(exons) <= end &
          as.character(S4Vectors::mcols(exons)[[name_col]]) %in% gene_names
  exons <- exons[keep]
  if (length(exons) == 0L)
    return(data.frame(gene = character(0), start = integer(0),
                      end = integer(0), strand = character(0),
                      stringsAsFactors = FALSE))
  data.frame(
    gene   = as.character(S4Vectors::mcols(exons)[[name_col]]),
    start  = pmax(start, GenomicRanges::start(exons)),
    end    = pmin(end,   GenomicRanges::end(exons)),
    strand = as.character(GenomicRanges::strand(exons)),
    stringsAsFactors = FALSE
  )
}

# Non-overlapping row assignment: greedy first-fit into rows so text labels don't collide.
.assign_gene_rows <- function(g_df) {
  n <- nrow(g_df)
  if (n == 0L) return(integer(0))
  ord <- order(g_df$start)
  rows <- integer(n)
  row_ends <- numeric(0)  # last end per row
  for (k in ord) {
    placed <- FALSE
    for (r in seq_along(row_ends)) {
      if (g_df$start[k] > row_ends[r]) {
        rows[k] <- r
        row_ends[r] <- g_df$end[k]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      rows[k] <- length(row_ends) + 1L
      row_ends <- c(row_ends, g_df$end[k])
    }
  }
  rows
}

.build_gene_track <- function(gene_df, exon_df, start, end, chrom, ax_fmt) {
  if (nrow(gene_df) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::geom_blank() +
        ggplot2::coord_cartesian(xlim = c(start, end)) +
        ggplot2::scale_x_continuous(labels = ax_fmt, expand = c(0, 0)) +
        ggplot2::labs(y = "Genes") +
        ggplot2::theme_classic(base_size = 10) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank())
    )
  }
  gene_df$row <- .assign_gene_rows(gene_df)

  # Attach row assignment to exons via gene name
  if (nrow(exon_df) > 0L) {
    exon_df$row <- gene_df$row[match(exon_df$gene, gene_df$gene)]
    exon_df <- exon_df[!is.na(exon_df$row), , drop = FALSE]
  }

  # Strand direction chevrons along the gene body (thin line)
  # ">" for + strand, "<" for - strand, placed every ~5% of region width
  span <- end - start
  step <- max(as.integer(span / 30), 1L)
  chev <- do.call(rbind, lapply(seq_len(nrow(gene_df)), function(i) {
    lo <- gene_df$start[i] + step
    hi <- gene_df$end[i]   - step
    if (hi <= lo) return(NULL)             # gene too short for chevrons
    xs <- seq(lo, hi, by = step)
    if (length(xs) == 0L) return(NULL)
    data.frame(
      x    = xs, y = gene_df$row[i],
      lab  = if (gene_df$strand[i] == "-") "<" else ">",
      stringsAsFactors = FALSE
    )
  }))

  gg <- ggplot2::ggplot() +
    # thin line for the whole gene body (intron backbone)
    ggplot2::geom_segment(
      data = gene_df,
      ggplot2::aes(x = .data$start, xend = .data$end,
                   y = .data$row,   yend = .data$row),
      color = "steelblue4", linewidth = 0.3
    )
  # strand chevrons on the intron line
  if (!is.null(chev) && nrow(chev) > 0L) {
    gg <- gg + ggplot2::geom_text(
      data = chev,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
      color = "steelblue4", size = 2.5
    )
  }
  # Thick exon rectangles
  if (nrow(exon_df) > 0L) {
    gg <- gg + ggplot2::geom_rect(
      data = exon_df,
      ggplot2::aes(xmin = .data$start, xmax = .data$end,
                   ymin = .data$row - 0.28, ymax = .data$row + 0.28),
      fill = "steelblue4", color = NA
    )
  }
  # Gene name labels above each gene row
  top_row <- max(gene_df$row)
  gg <- gg + ggplot2::geom_text(
      data = gene_df,
      ggplot2::aes(x = (.data$start + .data$end) / 2,
                   y = .data$row + 0.55,
                   label = .data$gene),
      size = 3, fontface = "italic"
    )

  gg +
    ggplot2::coord_cartesian(xlim = c(start, end),
                             ylim = c(0.3, top_row + 0.9)) +
    # Axis at TOP (IGV / UCSC ruler convention)
    ggplot2::scale_x_continuous(labels = ax_fmt, expand = c(0, 0),
                                position = "top") +
    # IGV-style stacked left-margin labels: bold `chr1` tag on top,
    # `Genes` axis label below.
    ggplot2::labs(y = "Genes", tag = chrom) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.text.y         = ggplot2::element_blank(),
      axis.ticks.y        = ggplot2::element_blank(),
      axis.title.y        = ggplot2::element_text(size = 10,
                                                  angle = 0, vjust = 0.5,
                                                  margin = ggplot2::margin(r = 4)),
      axis.title.x        = ggplot2::element_blank(),
      # bottom edge fully hidden (axis is on top now)
      axis.line.x.bottom  = ggplot2::element_blank(),
      axis.text.x.bottom  = ggplot2::element_blank(),
      axis.ticks.x.bottom = ggplot2::element_blank(),
      # top-left chrom tag (same size as the Genes axis label, unbold)
      plot.tag            = ggplot2::element_text(size = 10, hjust = 0),
      plot.tag.position   = c(0, 0.9)
    )
}

.build_snp_track <- function(snp_pos, snp_id, start, end, ax_fmt) {
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = snp_pos,
                        color = "#FF7917", linewidth = 0.6) +
    ggplot2::annotate("text", x = snp_pos, y = 0.5,
                      label = snp_id, hjust = -0.05,
                      size = 3, fontface = "bold") +
    ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(labels = ax_fmt, expand = c(0, 0)) +
    ggplot2::labs(y = "SNP") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x  = ggplot2::element_blank())
}

# Parabolic arcs (t in [0, 1]; y = 4*t*(1-t) scaled by loop span).
.arc_path <- function(x1, x2, n_pts = 60L, height_scale = 1) {
  t <- seq(0, 1, length.out = n_pts)
  span <- abs(x2 - x1)
  data.frame(
    x = x1 + t * (x2 - x1),
    y = 4 - (4 * t * (1 - t) * span * height_scale)
  )
}

.build_arc_panel <- function(arc_df, start, end, ax_fmt,
                              eff_lo, eff_hi, strong_color, weak_color,
                              panel_title, show_x = FALSE) {
  # Expand each arc into a path of points
  paths <- do.call(rbind, lapply(seq_len(nrow(arc_df)), function(i) {
    p <- .arc_path(arc_df$x1[i], arc_df$x2[i])
    p$loop_id <- arc_df$loop_id[i]
    p$value   <- arc_df$value[i]
    p
  }))
  y_max <- max(paths$y, na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  gg <- ggplot2::ggplot(paths, ggplot2::aes(x = .data$x, y = .data$y,
                                            group = .data$loop_id,
                                            color = .data$value)) +
    ggplot2::geom_path(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::scale_color_gradient(low = weak_color, high = strong_color,
                                  limits = c(eff_lo, eff_hi),
                                  oob = scales::squish,
                                  na.value = weak_color) +
    ggplot2::coord_cartesian(xlim = c(start, end),
                             ylim = c(0, y_max * 1.1)) +
    ggplot2::scale_x_continuous(labels = ax_fmt, expand = c(0, 0)) +
    ggplot2::labs(y = panel_title, color = "Mean strength") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y  = ggplot2::element_text(size = 9),
                   panel.grid = ggplot2::element_blank(),
                   # No X-axis line for the loops because it's at the top, above the genes
                   axis.line.x  = ggplot2::element_blank())
  return(gg)
}
