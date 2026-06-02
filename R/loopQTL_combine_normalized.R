#' Combine per-chromosome `loopQTL_normalize()` outputs into a full object
#'
#' Use this after running [loopQTL_normalize()] with disjoint `chroms` subsets
#' across multiple SLURM array tasks. Each per-chrom object has the same
#' `samples` / `loops` / `vcf_path` / `resolution`, only differing in which
#' rows of `phenotype` were filled and which chromosomes are listed in
#' `meta$chroms_normalized`. This function takes the per-chrom phenotype rows
#' from each input object and stitches them into a single phenotype matrix
#' covering all loops processed across the array.
#'
#' All input objects must share the same `samples`, `loops` (same length and
#' the same loop_ids), `resolution`, `vcf_path`, and `normalize_method`. The
#' processed chromosomes across inputs should be disjoint; if a loop's row is
#' filled in more than one input, the last input wins and a warning is issued.
#'
#' @param ... Two or more `loopQTL` objects produced by `loopQTL_normalize()`,
#'   or a single list of them.
#' @param paths Alternative input form: a character vector of `.rda` paths
#'   each containing an object named `obj` (or `x`); they're loaded and
#'   combined. Mutually exclusive with `...`.
#'
#' @return A single `loopQTL` object with `phenotype` rows filled for the
#'   union of all `chroms_normalized` across inputs. `meta$chroms_normalized`
#'   is the union; `meta$n_zero_var` is summed across inputs.
#'
#' @export
loopQTL_combine_normalized <- function(..., paths = NULL) {
  if (!is.null(paths)) {
    if (...length() > 0L)
      stop("Pass either `...` or `paths`, not both")
    objs <- lapply(paths, function(p) {
      e <- new.env()
      nms <- load(p, envir = e)
      # prefer `obj`, else `x`, else first loaded name
      if ("obj" %in% nms) e$obj
      else if ("x" %in% nms) e$x
      else get(nms[1L], envir = e)
    })
  } else {
    objs <- list(...)
    # allow caller to pass a single list
    if (length(objs) == 1L && is.list(objs[[1L]]) && !methods::is(objs[[1L]], "loopQTL"))
      objs <- objs[[1L]]
  }
  if (length(objs) < 1L)
    stop("No input objects provided")
  for (o in objs)
    if (!methods::is(o, "loopQTL"))
      stop("All inputs must be loopQTL objects")

  base <- objs[[1L]]
  base_loop_ids <- S4Vectors::mcols(base@loops)$loop_id
  n_loops <- length(base@loops)
  n_samples <- length(base@samples)

  # Validate compatibility + collect filled-row indicators
  for (k in seq_along(objs)) {
    o <- objs[[k]]
    if (!identical(o@samples, base@samples))
      stop("Object ", k, " has different `samples` from object 1")
    if (length(o@loops) != n_loops ||
        !identical(S4Vectors::mcols(o@loops)$loop_id, base_loop_ids))
      stop("Object ", k, " has a different consensus loop set from object 1")
    if (o@resolution != base@resolution)
      stop("Object ", k, " has different `resolution` from object 1")
    if (!identical(o@vcf_path, base@vcf_path))
      warning("Object ", k, " has a different `vcf_path` from object 1 (using object 1's)")
    if (length(o@phenotype) == 0L)
      stop("Object ", k, " has no phenotype filled; nothing to combine")
    if (!isTRUE(identical(o@meta$normalize_method, base@meta$normalize_method)))
      stop("Object ", k, " used normalize_method = ",
           deparse(o@meta$normalize_method), " but object 1 used ",
           deparse(base@meta$normalize_method))
  }

  # Build combined phenotype, last-writer-wins on overlap (with warning)
  pheno <- matrix(NA_real_, nrow = n_loops, ncol = n_samples,
                  dimnames = list(NULL, base@samples))
  fill_count <- integer(n_loops)
  for (o in objs) {
    filled_rows <- which(rowSums(!is.na(o@phenotype)) > 0L)
    if (length(filled_rows) == 0L) next
    overlap <- filled_rows[fill_count[filled_rows] > 0L]
    if (length(overlap) > 0L)
      warning("Overwriting ", length(overlap),
              " phenotype rows that were already filled by an earlier input")
    pheno[filled_rows, ] <- o@phenotype[filled_rows, , drop = FALSE]
    fill_count[filled_rows] <- fill_count[filled_rows] + 1L
  }

  base@phenotype <- pheno

  # ---- Merge cis ----
  # If any input flags user_provided, take that one as-is (genome-wide,
  # authoritative; all user_provided inputs should be identical anyway).
  # Otherwise sum per-chrom contributions across inputs.
  cis_sources <- vapply(objs,
                        function(o) as.character(o@meta$cis_source %||% NA_character_),
                        character(1L))
  user_inputs <- which(cis_sources == "user_provided")
  if (length(user_inputs) > 0L) {
    base@cis <- objs[[user_inputs[1L]]]@cis
    base@meta$cis_source <- "user_provided"
  } else if (all(cis_sources == "computed_from_hic", na.rm = TRUE) &&
             all(vapply(objs, function(o) length(o@cis) == n_samples, logical(1L)))) {
    cis_sum <- setNames(numeric(n_samples), base@samples)
    for (o in objs) cis_sum <- cis_sum + o@cis
    base@cis <- cis_sum
    base@meta$cis_source <- "computed_from_hic"
    base@meta$cis_chroms_in <- Reduce(union,
      lapply(objs, function(o) o@meta$cis_chroms_in %||% character(0)))
  } else {
    # Mixed / partial — leave whatever object 1 had and warn.
    warning("Could not merge cis cleanly across inputs (mixed sources or ",
            "incomplete cis vectors); object 1's cis preserved as-is")
  }

  # Merge other meta
  base@meta$chroms_normalized <- Reduce(union,
    lapply(objs, function(o) o@meta$chroms_normalized %||% character(0)))
  base@meta$combined_from <- length(objs)

  methods::validObject(base)
  base
}

# tiny helper
`%||%` <- function(a, b) if (is.null(a)) b else a
