#' Annotate SNP-loop QTL results with target genes, optional peak overlap,
#' and a Type classification.
#'
#' Adds three columns to `obj@results`:
#' \describe{
#'   \item{`gene`}{Comma-separated gene names whose TSS (± `TSS_window`) falls
#'     in the SNP's PARTNER anchor (anchor B if the SNP is in anchor A, and
#'     vice versa; both partners if `anchor == "both"`). Empty string if no
#'     gene is found.}
#'   \item{`peak`}{(only added when `peak` is supplied) `TRUE` if the SNP
#'     overlaps any peak from the supplied peak file, `FALSE` otherwise. When
#'     `remove_TSS = TRUE`, peaks that overlap any TSS region are dropped
#'     before this overlap check (so promoter-like peaks don't count).}
#'   \item{`Type`}{One of:
#'     \itemize{
#'       \item `"WithinTSS"` - SNP overlaps any gene's TSS (± `TSS_window`).
#'         Short-circuits all other annotation when `remove_TSS = TRUE`.
#'       \item `"rs-gene"` - SNP not in TSS, has a target gene at the partner anchor.
#'       \item `"rs-<peak_annotation>"` - SNP not in TSS, no gene at partner,
#'         but partner anchor overlaps a peak (only when `peak` is supplied).
#'         e.g. `"rs-Enhancer"`.
#'       \item `"rs-O"` - SNP not in TSS, no gene AND no peak at the partner anchor.
#'       \item `"rs(<peak_annotation>)-gene"` - additionally overlaps a peak
#'         (only when `peak` is supplied), has a target gene. e.g.
#'         `"rs(Enhancer)-gene"`.
#'       \item `"rs(<peak_annotation>)-<peak_annotation>"` - SNP overlaps peak,
#'         partner overlaps peak, no gene at partner.
#'         e.g. `"rs(Enhancer)-Enhancer"`.
#'       \item `"rs(<peak_annotation>)-O"` - overlaps a peak but no target gene
#'         and no peak at partner.
#'     }
#'     Partner-side priority: `gene` > peak > `O`.
#'   }
#' }
#'
#' @param obj A [loopQTL-class] with `obj@results` populated by `runQTL()`.
#' @param gtf Path to a GTF/GFF file with `gene` entries. TSS is taken as
#'   the gene's `start` on + strand, `end` on - strand. Read via
#'   `rtracklayer::import()`.
#' @param peak Optional path to a peak file (BED-like; anything
#'   `rtracklayer::import()` can read). If supplied, an extra `peak` column
#'   is added and the `Type` labels gain the `"rs(<peak_annotation>)-..."` form.
#' @param TSS_window Bp flank around each gene's TSS. Default 2000. Used both
#'   to define "TSS regions" (for `WithinTSS` classification and peak
#'   filtering when `remove_TSS = TRUE`) and to count a gene as a target
#'   when its TSS region falls within `anchor_window` of the partner anchor.
#' @param anchor_window Bp flank added around each anchor when looking for
#'   overlapping TSSs (and peaks). Default 5000.
#' @param peak_annotation Character label substituted into the `Type` values
#'   when `peak` is supplied. Default `"Peak"` -> `"rs(Peak)-gene"`; common
#'   alternatives are `"Enhancer"` or `"H3K27ac"`.
#' @param remove_TSS Logical; default `TRUE`. When `TRUE`:
#'   \enumerate{
#'     \item SNPs overlapping any TSS region get `Type = "WithinTSS"`,
#'       short-circuiting their gene/peak annotation.
#'     \item If `peak` is supplied, peaks that overlap any TSS region are
#'       dropped before the SNP-in-peak overlap check (so the `peak` column
#'       reflects only TSS-distal peaks).
#'   }
#'   When `FALSE`, no TSS-based filtering is done on either side.
#'
#' @return The input `obj` with `obj@results` augmented by `gene`, `peak`
#'   (when applicable), and `Type` columns.
#'
#' @export
loopQTLannotate <- function(obj,
                             gtf,
                             peak            = NULL,
                             TSS_window      = 2000L,
                             anchor_window   = 5000L,
                             peak_annotation = "Peak",
                             remove_TSS      = TRUE) {
  stopifnot(
    methods::is(obj, "loopQTL"),
    is.character(gtf), length(gtf) == 1L, file.exists(gtf),
    is.numeric(TSS_window),    length(TSS_window)    == 1L, TSS_window    >= 0,
    is.numeric(anchor_window), length(anchor_window) == 1L, anchor_window >= 0,
    is.character(peak_annotation), length(peak_annotation) == 1L,
    is.logical(remove_TSS), length(remove_TSS) == 1L
  )
  if (!is.null(peak)) {
    stopifnot(is.character(peak), length(peak) == 1L, file.exists(peak))
  }
  if (nrow(obj@results) == 0L)
    stop("obj@results is empty. Run runQTL(obj) first.")

  TSS_window    <- as.integer(TSS_window)
  anchor_window <- as.integer(anchor_window)

  # 1. Read GTF, extract gene entries, derive TSS
  message("[", format(Sys.time()), "] reading GTF: ", gtf)
  gtf_gr <- rtracklayer::import(gtf)
  type_col <- if ("type" %in% names(S4Vectors::mcols(gtf_gr))) "type" else
              if ("feature" %in% names(S4Vectors::mcols(gtf_gr))) "feature" else
              NULL
  if (is.null(type_col))
    stop("GTF has neither `type` nor `feature` column; cannot find gene entries")
  gene_gr <- gtf_gr[S4Vectors::mcols(gtf_gr)[[type_col]] == "gene"]
  if (length(gene_gr) == 0L)
    stop("No `gene` entries in GTF")

  gene_names <- if ("gene_name" %in% names(S4Vectors::mcols(gene_gr))) {
    as.character(gene_gr$gene_name)
  } else if ("gene_id" %in% names(S4Vectors::mcols(gene_gr))) {
    as.character(gene_gr$gene_id)
  } else {
    stop("GTF gene entries lack both `gene_name` and `gene_id`")
  }
  # Some GTFs have NA gene_name for novel/predicted; fall back to gene_id.
  if ("gene_id" %in% names(S4Vectors::mcols(gene_gr))) {
    na_idx <- which(is.na(gene_names) | gene_names == "")
    if (length(na_idx) > 0L)
      gene_names[na_idx] <- as.character(gene_gr$gene_id[na_idx])
  }

  strand_chr <- as.character(GenomicRanges::strand(gene_gr))
  tss_pos <- ifelse(strand_chr == "-",
                    GenomicRanges::end(gene_gr),
                    GenomicRanges::start(gene_gr))
  tss_regions <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gene_gr),
    ranges   = IRanges::IRanges(start = pmax(1L, tss_pos - TSS_window),
                                end   = tss_pos + TSS_window)
  )
  S4Vectors::mcols(tss_regions)$gene <- gene_names
  message("[", format(Sys.time()), "] ", length(tss_regions),
          " gene TSS regions (±", TSS_window, " bp)")

  # 2. Chrom-naming normalizer for THIS GTF  # GTF may use "chr1" (GENCODE) or "1" (Ensembl); loop anchors may use
  # either depending on the BEDPE caller (Mustache uses "1"). Build a
  # function that converts any chrom string to whatever the GTF uses.
  gtf_chroms <- unique(as.character(GenomicRanges::seqnames(tss_regions)))
  gtf_uses_chr <- any(grepl("^chr", gtf_chroms))
  normalize_chrom <- function(x) {
    x <- as.character(x)
    if (gtf_uses_chr) ifelse(grepl("^chr", x), x, paste0("chr", x))
    else              sub("^chr", "", x)
  }

  # 3. Optional: read peaks, filter by TSS if remove_TSS
  peaks_gr <- NULL
  if (!is.null(peak)) {
    message("[", format(Sys.time()), "] reading peaks: ", peak)
    peaks_gr <- rtracklayer::import(peak)
    # Normalize peak chroms to GTF style so overlap actually works.
    peak_chroms_norm <- normalize_chrom(GenomicRanges::seqnames(peaks_gr))
    peaks_gr <- GenomicRanges::GRanges(
      seqnames = peak_chroms_norm,
      ranges   = IRanges::ranges(peaks_gr)
    )
    if (isTRUE(remove_TSS)) {
      hits <- GenomicRanges::findOverlaps(peaks_gr, tss_regions)
      peaks_in_tss <- unique(S4Vectors::queryHits(hits))
      n_dropped <- length(peaks_in_tss)
      if (n_dropped > 0L) {
        peaks_gr <- peaks_gr[-peaks_in_tss]
        message("[", format(Sys.time()), "] dropped ", n_dropped,
                " peaks overlapping TSS regions; ", length(peaks_gr),
                " peaks retained for SNP overlap.")
      }
    }
  }

  # 4. Per-loop gene lists at anchor A and anchor B  # Critical: normalize anchor seqnames to the GTF's convention BEFORE
  # findOverlaps. Without this, anchor seqnames like "1" silently fail to
  # match TSS seqnames like "chr1", producing zero gene assignments. (This
  # was the cause of the original "no rs-gene at all" bug.)
  loops_gi <- obj@loops
  loop_ids_all <- S4Vectors::mcols(loops_gi)$loop_id
  anchorA <- InteractionSet::anchors(loops_gi, type = "first")
  anchorB <- InteractionSet::anchors(loops_gi, type = "second")

  pad_anchor <- function(gr) {
    GenomicRanges::GRanges(
      seqnames = normalize_chrom(GenomicRanges::seqnames(gr)),
      ranges   = IRanges::IRanges(
        start = pmax(1L, GenomicRanges::start(gr) - anchor_window),
        end   = GenomicRanges::end(gr) + anchor_window
      )
    )
  }
  padded_A <- pad_anchor(anchorA)
  padded_B <- pad_anchor(anchorB)

  # Anchor-side gene assignments (length = n_loops; one character vector per loop)
  build_gene_lists <- function(padded) {
    hits <- GenomicRanges::findOverlaps(padded, tss_regions)
    out <- vector("list", length(padded))
    if (length(hits) == 0L) return(out)
    by_loop <- split(S4Vectors::subjectHits(hits),
                     S4Vectors::queryHits(hits))
    for (lp in names(by_loop)) {
      out[[as.integer(lp)]] <-
        unique(gene_names[by_loop[[lp]]])
    }
    out
  }
  genes_at_A <- build_gene_lists(padded_A)
  genes_at_B <- build_gene_lists(padded_B)

  # Per-anchor peak overlap (logical vector, length = n_loops). NULL when no
  # peak file was supplied -- then partner-side never gets a peak label.
  peak_at_A <- if (!is.null(peaks_gr)) {
    GenomicRanges::countOverlaps(padded_A, peaks_gr) > 0L
  } else NULL
  peak_at_B <- if (!is.null(peaks_gr)) {
    GenomicRanges::countOverlaps(padded_B, peaks_gr) > 0L
  } else NULL

  # 5. Per-row annotation in obj@results
  res <- obj@results
  loop_idx <- match(res$loop_id, loop_ids_all)

  # SNP GRanges (1 bp per SNP), with chrom normalized to GTF convention.
  snp_gr <- GenomicRanges::GRanges(
    seqnames = normalize_chrom(res$chrom),
    ranges   = IRanges::IRanges(start = res$pos, end = res$pos)
  )

  snp_in_tss <- if (isTRUE(remove_TSS)) {
    GenomicRanges::countOverlaps(snp_gr, tss_regions) > 0L
  } else {
    rep(FALSE, nrow(res))
  }
  snp_in_peak <- if (!is.null(peaks_gr)) {
    GenomicRanges::countOverlaps(snp_gr, peaks_gr) > 0L
  } else {
    rep(NA, nrow(res))
  }

  # Resolve partner-anchor genes per row
  partner_genes <- vector("list", nrow(res))
  partner_has_peak <- rep(FALSE, nrow(res))
  for (i in seq_len(nrow(res))) {
    li <- loop_idx[i]
    anc <- res$anchor[i]
    g <- character(0)
    if (anc == "A")         g <- genes_at_B[[li]]
    else if (anc == "B")    g <- genes_at_A[[li]]
    else if (anc == "both") g <- unique(c(genes_at_A[[li]], genes_at_B[[li]]))
    partner_genes[[i]] <- if (is.null(g)) character(0) else g
    if (!is.null(peak_at_A)) {
      partner_has_peak[i] <-
        if (anc == "A")         peak_at_B[li]
        else if (anc == "B")    peak_at_A[li]
        else                    (peak_at_A[li] || peak_at_B[li])
    }
  }

  gene_col <- vapply(partner_genes,
                     function(g) paste(g, collapse = ","),
                     character(1L))
  has_gene <- nchar(gene_col) > 0L

  # Build Type column. Partner-side priority: gene > peak > O.
  rs_prefix <- if (!is.null(peaks_gr)) {
    ifelse(snp_in_peak, paste0("rs(", peak_annotation, ")"), "rs")
  } else {
    rep("rs", nrow(res))
  }
  partner_suffix <- ifelse(
    has_gene, "gene",
    ifelse(partner_has_peak, peak_annotation, "O")
  )
  type_col <- paste0(rs_prefix, "-", partner_suffix)
  type_col[snp_in_tss] <- "WithinTSS"

  res$gene <- gene_col
  if (!is.null(peaks_gr)) res$peak <- snp_in_peak
  res$Type <- type_col

  obj@results <- res
  obj@meta$annotate <- list(
    gtf            = gtf,
    peak           = peak,
    TSS_window     = TSS_window,
    anchor_window  = anchor_window,
    peak_annotation = peak_annotation,
    remove_TSS     = remove_TSS,
    n_genes_used   = length(unique(gene_names))
  )
  obj
}
