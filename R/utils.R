# Internal helpers. Not exported.

# Bin a 1-based BEDPE-style anchor (start, end) to its analysis-bin start (0-based),
# matching strawr's binning. straw returns bin starts as the lower coordinate of
# the bin in 0-based half-open BP space, so bin = floor(midpoint / resolution) * resolution.
.bin_anchor <- function(start, end, resolution) {
  mid <- (as.numeric(start) + as.numeric(end)) / 2
  as.integer(floor(mid / resolution) * resolution)
}

# Detect strawr's chromosome naming style ("1" vs "chr1") from the .hic header.
# Returns a function that maps a user-supplied chrom string (either style) to the
# style this .hic expects.
.make_chrom_normalizer <- function(hic_chroms) {
  uses_chr_prefix <- any(grepl("^chr", hic_chroms))
  function(x) {
    if (uses_chr_prefix) {
      ifelse(grepl("^chr", x), x, paste0("chr", x))
    } else {
      sub("^chr", "", x)
    }
  }
}

# RINT with Blom offset, NA-aware. Operates on a single vector (one loop's
# values across samples).
.rint <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  n <- sum(!is.na(r))
  if (n < 2L) return(rep(NA_real_, length(x)))
  qnorm((r - 0.5) / n)
}

# Apply .rint() row-wise to a (loops x samples) matrix. Rows whose values are
# all-NA, fewer than 2 non-NA, or have zero variance are returned as all-NA
# (RINT on a constant produces NaN/Inf, which would poison downstream
# regressions). dimnames are preserved.
.rint_matrix <- function(mat) {
  if (length(mat) == 0L) return(mat)
  out <- matrix(NA_real_,
                nrow = nrow(mat), ncol = ncol(mat),
                dimnames = dimnames(mat))
  row_sd <- apply(mat, 1L, stats::sd, na.rm = TRUE)
  ok <- is.finite(row_sd) & row_sd > 0
  if (any(ok)) {
    good <- which(ok)
    out[good, ] <- t(apply(mat[good, , drop = FALSE], 1L, .rint))
  }
  out
}

# Read a single BEDPE-like loop file. Tolerant of "#chr" headers and of
# files with or without a `loop_id` column.
.read_loops_bedpe <- function(path) {
  if (!file.exists(path))
    stop("Loop file not found: ", path)
  hdr <- readLines(path, n = 1L, warn = FALSE)
  has_header <- grepl("^#?chr", hdr, ignore.case = TRUE) &&
                grepl("\\bstart", hdr, ignore.case = TRUE)
  df <- utils::read.table(
    path, header = has_header, sep = "\t",
    comment.char = if (has_header) "" else "#",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!has_header) {
    # assume positional: chrA startA endA chrB startB endB [loop_id]
    if (ncol(df) < 6L)
      stop("Loop file ", path, " must have at least 6 columns")
    base_names <- c("chrA","startA","endA","chrB","startB","endB")
    if (ncol(df) == 6L) {
      names(df) <- base_names
    } else {
      names(df) <- c(base_names, paste0("V", seq_len(ncol(df) - 6L) + 6L))
    }
  } else {
    nm <- names(df)
    if (nm[1] %in% c("#chr1", "#chrA")) nm[1] <- "chrA"
    # canonicalize common variants
    nm <- sub("^chr1$",  "chrA",   nm)
    nm <- sub("^x1$",    "startA", nm, ignore.case = TRUE)
    nm <- sub("^y1$",    "endA",   nm, ignore.case = TRUE)
    nm <- sub("^chr2$",  "chrB",   nm)
    nm <- sub("^x2$",    "startB", nm, ignore.case = TRUE)
    nm <- sub("^y2$",    "endB",   nm, ignore.case = TRUE)
    names(df) <- nm
  }
  need <- c("chrA","startA","endA","chrB","startB","endB")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0L)
    stop("Loop file ", path, " missing required columns: ",
         paste(miss, collapse = ", "))
  df$chrA <- as.character(df$chrA)
  df$chrB <- as.character(df$chrB)
  df$startA <- as.integer(df$startA)
  df$endA   <- as.integer(df$endA)
  df$startB <- as.integer(df$startB)
  df$endB   <- as.integer(df$endB)
  df
}
