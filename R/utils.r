#' Calculate Delta betas between two groups
#'
#' @param betas array of betas values
#' @param design design matrix
#' @param cmtx contrast matrix
#' @param contrast_name column
#'
#' @importFrom MatrixGenerics rowMeans
#'
#' @return vector
get_delta_betas <- function(betas, design, cmtx, contrast_name) {

  score <- round(as.vector(design %*% cmtx[, contrast_name]), 6)

  group1 <- score > 0
  group0 <- score < 0

  group1_means <- MatrixGenerics::rowMeans(
    betas[, group1, drop = FALSE],
    na.rm = TRUE
  )
  group0_means <- MatrixGenerics::rowMeans(
    betas[, group0, drop = FALSE],
    na.rm = TRUE
  )

  (group1_means - group0_means) * 100
}


##' Calculate Delta betas between two groups
##'
##' @param betas array of betas values
##' @param design design matrix
##' @param factor_prefix prefix for factor columns
##' @param level_col column for the level of interest
##'
##' @importFrom MatrixGenerics rowMeans
##'
##' @return vector
# get_delta_betas <- function(betas, design, factor_prefix, level_col) {
#
#  factor_cols <- grep(
#    paste0("^", factor_prefix),
#     colnames(design),
#     value = TRUE
#   )

#   group1 <- design[, level_col] == 1
#   group0 <- rowSums(design[, factor_cols, drop = FALSE]) == 0

#   group1_means <- MatrixGenerics::rowMeans(
#     betas[, group1, drop = FALSE],
#     na.rm = TRUE
#   )

#   group0_means <- MatrixGenerics::rowMeans(
#     betas[, group0, drop = FALSE],
#     na.rm = TRUE
#   )

#   (group1_means - group0_means) * 100
# }

#' Calculate mvalues
#'
#' @param betas array of betas values
#'
#' @return matrix
beta2m <- function(betas) {
  log2(betas / (1 - betas))
}

#' Calculate mvalues
#'
#' @param mvalues array of M values
#'
#' @return matrix
#'
m2beta <- function(mvalues) {
  (2^mvalues) / (1 + 2^mvalues)
}

#' Winsorize beta values independently for each probe
#'
#' Caps values in each probe row at row-specific quantiles.
#'
#' @param betas Numeric matrix of beta values, with probes in rows and samples
#'   in columns.
#' @param probs Numeric vector of length two giving the lower and upper
#'   quantiles used as caps (default \code{c(0.01, 0.99)}).
#' @return A matrix with the same dimensions and dimnames as \code{betas};
#'   missing values are preserved.
#' @noRd
winsorize_by_probe <- function(betas, probs = c(0.01, 0.99)) {
  if (!is.matrix(betas) || !is.numeric(betas)) {
    stop("`betas` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      any(!is.finite(probs)) || probs[1] < 0 || probs[2] > 1 ||
      probs[1] > probs[2]) {
    stop(
      "`probs` must contain two ordered probabilities between 0 and 1.",
      call. = FALSE
    )
  }

  for (i in seq_len(nrow(betas))) {
    observed <- !is.na(betas[i, ])
    if (!any(observed)) {
      next
    }

    limits <- stats::quantile(
      betas[i, observed],
      probs = probs,
      names = FALSE
    )
    betas[i, observed] <- pmax(
      limits[1],
      pmin(betas[i, observed], limits[2])
    )
  }

  betas
}


#' Filter probes from list of probes
#'
#' @param filters file containing probe list
#' @param plateform plateform
#' @return vector
#'
#' @importFrom RCurl url.exists
#' @importFrom readr read_tsv
cpg_excl <- function(
  filters = NULL,
  plateform = "IlluminaHumanMethylationEPIC",
  snp = TRUE,
  cross = TRUE,
  xy = TRUE
) {

  if (is.null(filters)) {
    # from list
    return(get_default_probeList(plateform, snp, cross, xy))
  }
  # from files
  probes <- c()
  for (file in filters) {
    if (!file.exists(file) && !RCurl::url.exists(file)) {
      stop(paste0(file, " not found"))
    }
    probes <- unique(c(probes, readr::read_tsv(file)[[1]]))
  }
  probes
}

#' Generic function to extract column data from an object
#'
#' @param x An object RGset or sdfs list.
#' @param outfile file output name
#'
#' @export
#'
setGeneric("to_geo_submission",
  function(
    x,
    outfile,
    digits = 4
  ) {
    standardGeneric("to_geo_submission")
  }
)


#' generate files to submit data to Geo database
#'
#' @param RGset an RGset object
#' @param outfile file output name
#'
#' @importFrom readr write_tsv
#' @importFrom dplyr as_tibble
#'
#' @export
setMethod(
  "to_geo_submission",
  signature(x = "RGChannelSet", outfile = "character", digits = "numeric"),
  definition = function(x, outfile = "rawdata2Geo.tsv", digits = 4) {

    if (!requireNamespace("minfi", quietly = TRUE)) {
      stop("Package 'minfi' is required for this function.",
        "Please install it.",
        call. = FALSE
      )
    }
    mset <- minfi::preprocessRaw(x)
    meth <- minfi::getMeth(mset)
    colnames(meth) <-
      paste(minfi::pData(x)$Basename, "Methylated signal")
    unmeth <- minfi::getUnmeth(mset)
    colnames(unmeth) <-
      paste(minfi::pData(x)$Basename, "Unmethylated signal")
    pval <- minfi::detectionP(x)
    colnames(pval) <- paste(minfi::pData(x)$Basename, "Detection Pval")

    msi <- cbind(unmeth, meth, pval)
    n <- ncol(pval) + 1
    m <- ncol(pval) * 2 + 1
    cols <- c(1, n, m) + rep(0:(n - 2), each = 3)
    msi <- msi[, cols]
    readr::write_tsv(dplyr::as_tibble(msi), quote = FALSE, file = outfile)

  }
)

#' generate files to submit data to Geo database
#'
#' @param Betas a Betas object
#' @param outfile file output name
#' @param digits number of decimal places to round to
#'
#' importFrom assertthat asserthat
#' importFrom readr write_tsv
#'
#' @export
setMethod(
  "to_geo_submission",
  signature(x = "list", outfile = "character", digits = "numeric"),
  definition = function(x, outfile = "betas2Geo.tsv", digits = 4) {

    if (is.null(names(x)) || any(names(x) == "")) {
      stop("`sdfs` must be a named list")
    }

    # --- betas and detection p-values, extract from sdf ---
    betas_list <- lapply(x, getBetas)
    pvals_list <- lapply(x,
      function(sdf) sesame::pOOBAH(sdf, return.pval = TRUE)
    )

    # same set of probes for everyone
    # (union, not intersection: we keep everything)
    all_probes <- Reduce(union, lapply(betas_list, names))

    betas_matrix <- sapply(betas_list, function(x) x[all_probes])
    pvals_matrix <- sapply(pvals_list, function(x) x[all_probes])
    rownames(betas_matrix) <- all_probes
    rownames(pvals_matrix) <- all_probes

    # security : same columns, same order
    pvals_matrix <- pvals_matrix[, colnames(betas_matrix), drop = FALSE]

    # --- construction dof combined format ---
    combined <- data.frame(ID_REF = all_probes, check.names = FALSE)
    for (s in colnames(betas_matrix)) {
      combined[[s]] <- round(betas_matrix[, s], digits)
      combined[[paste0(s, ".Detection Pval")]] <-
        round(pvals_matrix[, s], digits)
    }

    write.table(
      combined,
      file = outfile,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    invisible(combined)
  }
)

#' h_mean
#' Calculate the harmonic mean of a numeric vector
#' @param x A numeric vector
#' @param na.rm Logical, whether to remove NA values
#'   before calculation (default: TRUE)
#' @return The harmonic mean of the input vector
h_mean <- function(x, na.rm = TRUE) {
  # Suppress NA si asked
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)

  # To avoid division by zero
  if (any(x == 0)) return(0)

  # Calcul classique
  (1 / mean(1 / x))
}


#' Regroupe les DMRs qui partagent au moins une sonde
#'
#' @param edges data.frame avec les colonnes `ID` (DMR) et `Probe_ID`
#'   (une ligne par couple DMR / sonde).
#' @return Vecteur nomme (noms = ID de DMR, valeurs = numero de cluster),
#'   ou NULL si aucune arete valide.
#' @noRd
dmr_clusters <- function(edges) {
  edges <- edges |>
    dplyr::transmute(
      ID = trimws(as.character(.data$ID)),
      Probe_ID = trimws(as.character(.data$Probe_ID))
    ) |>
    dplyr::filter(
      !is.na(.data$ID), .data$ID != "",
      !is.na(.data$Probe_ID), .data$Probe_ID != ""
    ) |>
    dplyr::distinct()

  if (nrow(edges) == 0) {
    return(NULL)
  }

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  igraph::components(g)$membership
}



#' Valide et compile une expression `combine`
#'
#' Grammaire : noms d'outils, `&`, `|`, `!`, `( )`, `n_of(k, a, b, ...)`
#' (au moins k outils parmi ceux mentionnes) et `only(a, b, ...)` (exactement
#' ces outils parmi ceux mentionnes). Tout le reste est refuse.
#'
#' @param expr chaine de caracteres.
#' @param tools noms d'outils valides (univers de `only()`).
#' @return Une fonction `f(pres)` qui prend une matrice logique
#'   regions x outils et renvoie un vecteur logique. L'attribut `tools` de la
#'   fonction contient les outils references dans l'expression.
#' @noRd
compile_combine <- function(expr, tools) {
  if (!is.character(expr) || length(expr) != 1 || is.na(expr)) {
    stop("`combine` must be a single character string.", call. = FALSE)
  }

  parsed <- tryCatch(
    str2lang(expr),
    error = function(e) {
      stop("Invalid `combine` expression: ", conditionMessage(e), call. = FALSE)
    }
  )

  bad_tool <- function(nm) {
    stop(
      "Unknown tool '", nm, "' in `combine`. Available: ",
      paste(tools, collapse = ", "), call. = FALSE
    )
  }

  collect_tools <- function(node) {
    if (is.symbol(node)) {
      return(as.character(node))
    }
    if (is.call(node)) {
      return(unlist(lapply(as.list(node)[-1], collect_tools),
                    use.names = FALSE))
    }
    character()
  }
  referenced_tools <- unique(collect_tools(parsed))

  walk <- function(node) {
    if (is.symbol(node)) {
      nm <- as.character(node)
      if (!nm %in% tools) bad_tool(nm)
      return(node)
    }

    if (!is.call(node) || !is.symbol(node[[1]])) {
      stop("Unsupported element in `combine`: ", deparse(node), call. = FALSE)
    }

    fn <- as.character(node[[1]])
    args <- as.list(node)[-1]

    if (fn %in% c("&", "|", "!", "(")) {
      node[-1] <- lapply(args, walk)
      return(node)
    }

    if (fn == "n_of") {
      if (
        length(args) < 2 || !is.numeric(args[[1]]) || length(args[[1]]) != 1
      ) {
        stop("n_of() expects a number followed by tool names: n_of(2, a, b, c).",
             call. = FALSE)
      }
      node[-c(1, 2)] <- lapply(args[-1], walk)
      return(node)
    }

    if (fn == "only") {
      if (length(args) < 1 || !all(vapply(args, is.symbol, logical(1)))) {
        stop("only() expects tool names: only(a, b).", call. = FALSE)
      }
      chosen <- vapply(args, as.character, character(1))
      invisible(lapply(args, walk))
      others <- setdiff(unique(referenced_tools), chosen)
      terms <- c(
        lapply(chosen, as.symbol),
        lapply(others, function(o) call("!", as.symbol(o)))
      )
      return(Reduce(function(a, b) call("&", a, b), terms))
    }

    stop("Function '", fn, "' is not allowed in `combine`.", call. = FALSE)
  }

  ast <- walk(parsed)

  predicate <- function(pres) {
    env <- new.env(parent = emptyenv())
    for (op in c("&", "|", "!", "(")) {
      env[[op]] <- get(op, envir = baseenv())
    }
    env$n_of <- function(k, ...) rowSums(cbind(...)) >= k
    for (tl in unique(referenced_tools)) {
      env[[tl]] <- unname(pres[, tl])
    }
    res <- eval(ast, env)
    res & !is.na(res)
  }
  attr(predicate, "tools") <- unique(referenced_tools)
  predicate
}
