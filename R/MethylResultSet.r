#' Define a Class Union
#'
#' Defines a class union accepting either a list or NULL. Used for optional
#' slots in S4 classes.
#'
#' @keywords classes
#' @rdname setClassUnion
setClassUnion("List_OR_NULL", c("list", "NULL"))

#' MethylResultSet Class
#'
#' An S4 class for storing and working with results from Illumina methylation
#' array analysis. Produced by fitting a linear model with limma on M-values.
#'
#' @slot manifest A \code{DFrame} containing the array manifest.
#' @slot dmps A \code{List_OR_NULL} where each element is a data frame of
#'   differential methylation positions (DMPs) for one model contrast.
#' @slot lambda A \code{numeric} vector of genomic inflation factors.
#' @slot metadata A \code{list} containing analysis metadata (model formula,
#'   intercept group, fitting method, and any metadata inherited from the input
#'   \code{SummarizedExperiment}).
#'
#' @docType class
#' @name MethylResultSet
#' @aliases MethylResultSet-class
#'
#' @rdname MethylResultSet
#'
#' @importFrom S4Vectors DataFrame
#'
#' @export
#'
MethylResultSet <- setClass(
  "MethylResultSet",
  representation(
    manifest  = "DFrame",
    dmps      = "List_OR_NULL",
    dmrs      = "List_OR_NULL",
    lambda    = "numeric",
    metadata  = "list"
  ),
  prototype(
    manifest  = S4Vectors::DataFrame(),
    dmps      = list(),
    dmrs      = list(),
    lambda    = numeric(0),
    metadata  = list()
  )
)

#' Constructor for MethylResultSet
#'
#' Fits a limma linear model on M-values from a \code{SummarizedExperiment}
#' and returns a \code{MethylResultSet} object containing per-probe statistics
#' for every model contrast.
#'
#' @param se A \code{SummarizedExperiment} (or subclass) produced by
#'   \code{methylkey}. Must contain M-values accessible via \code{get_mvals()}
#'   and sample metadata via \code{colData()}.
#' @param model A model formula as a character string, e.g. \code{"~group"} or
#'   \code{"~group+age"}. The first variable after \code{~} is used as the
#'   grouping factor.
#' @param intercept The reference level for the grouping factor (character).
#' @param method Fitting method passed to \code{limma::lmFit}. One of
#'   \code{"ls"} (default) or \code{"robust"}.
#'
#' @return An object of class \code{\link{MethylResultSet}}.
#'
#' @seealso \code{\link{MethylResultSet-class}}, \code{\link{getDMPs}},
#'   \code{\link{get_results}}
#'
#' @importFrom limma makeContrasts contrasts.fit eBayes
#' @importFrom stats as.formula model.matrix relevel qchisq
#' @importFrom dplyr mutate select arrange left_join
#' @importFrom tibble rownames_to_column
#'
#' @rdname MethylResultSet
#' @export
MethylResultSet <- function(
    se,
    model,
    intercept,
    method  = "ls") {

  pdata <- data.frame(colData(se))
  mval <- get_mvals(se)
  model <- tolower(model)

  # Extract the primary grouping variable robustly
  grp_g <- all.vars(stats::as.formula(model))[1]
  pdata[, grp_g] <- relevel(as.factor(unlist(pdata[, grp_g])), intercept)

  formula1 <- stats::as.formula(model)
  design   <- stats::model.matrix(formula1, data  = pdata)
  colnames(design) <- make.names(colnames(design))
  cmtx <- limma::makeContrasts(
    contrasts = colnames(design),
    levels  = colnames(design)
  )

  # Fit linear model
  fit <- limma::lmFit(mval, design, pdata, ndups  = 1, method  = method)
  rownames(cmtx) <- colnames(fit)
  fit_contrasts <- limma::contrasts.fit(fit, cmtx)
  eb <- limma::eBayes(fit_contrasts)

  # Genomic inflation factor (lambda) per contrast
  chisq  <- stats::qchisq(1 - eb$p.value, 1)
  lambda <- apply(chisq, 2, function(x) {
    stats::median(x, na.rm = TRUE) / stats::qchisq(0.5, 1)
  })

  # R-squared (goodness of fit) per probe
  sst <- rowSums(mval^2)
  ssr <- sst - fit$df.residual * fit$sigma^2
  rsq <- ssr / sst

  dmps <- lapply(colnames(fit), function(x) {
    top_tables(eb, x, rsq)
  })
  names(dmps) <- colnames(fit)

  for (x in names(dmps)) {
    #dmps[[x]]$deltabetas <- get_delta_betas(
    #  m2beta(mval[rownames(dmps[[x]]), ]),
    #  design[, x]
    #)

    dmps[[x]]$deltabetas <- get_delta_betas(
      betas = m2beta(mval[rownames(dmps[[x]]), ]),
      design = design,
      factor_prefix = grp_g,
      level_col = x
    )

    dmps[[x]] <- dmps[[x]] |>
      dplyr::mutate(status  = ifelse(.data$deltabetas > 0, "hyper", "hypo")) |>
      tibble::rownames_to_column("Probe_ID") |>
      dplyr::select(-"logFC") |>
      dplyr::arrange(.data$Probe_ID)
  }

  metadata <- S4Vectors::metadata(se)
  metadata$model <- model
  metadata$intercept <- intercept
  metadata$method   <- method


  new("MethylResultSet",
    manifest  = get_annotated_manifest(se),
    dmps      = dmps,
    lambda    = lambda,
    metadata  = metadata
  )
}

#' Retrieve Differentially Methylated Positions (DMPs)
#'
#' Returns the DMP table for a given contrast, joined with probe annotations
#' from the array manifest.
#'
#' @param x A \code{MethylResultSet} object.
#' @param group Integer index or character name of the contrast to retrieve
#'   (default \code{1}).
#'
#' @return A tibble with one row per probe, containing limma statistics and
#'   manifest annotations.
#'
#' @export
setGeneric("get_dmps", function(x, group  = 1) {
  standardGeneric("get_dmps")
})

#' @importFrom dplyr left_join rename mutate select
#' @describeIn get_dmps Method for MethylResultSet
setMethod("get_dmps", "MethylResultSet",
  function(x, group  = 1) {
    if (is.numeric(group) && (group < 1 || group > length(x@dmps))) {
      warning(sprintf(
        "group index %d is out of bounds (MethylResultSet has %d contrast(s)).",
        "Returning NULL.", group, length(x@dmps)
      ))
      return(NULL)
    }
    dplyr::left_join(
      x@dmps[[group]],
      dplyr::as_tibble(x@manifest),
      by  = "Probe_ID"
    ) |>
      dplyr::rename(
        chr = "CpG_chrm",
        strand = "probe_strand"
      ) |>
      dplyr::mutate(pos = .data$CpG_beg + 1) |>
      dplyr::select(-c("CpG_beg", "CpG_end"))
  }
)

#' Retrieve DMP Ranges as a GRanges Object
#'
#' Returns significant DMPs for a given contrast as a \code{GRanges} object,
#' filtered by adjusted p-value.
#'
#' @param x A \code{MethylResultSet} object.
#' @param index Integer index or character name of the contrast
#' (default \code{1}).
#' @param q Adjusted p-value threshold (default \code{0.05}).
#'
#' @return A \code{GRanges} object with one range per significant probe.
#'
#' @export
setGeneric("get_dmps_ranges",
  function(
    x,
    index = NULL,
    q = 0.05,
    keep_extra_columns = TRUE,
    as_list = FALSE
  ) {
    standardGeneric("get_dmps_ranges")
  }
)

#' @describeIn get_dmps_ranges Method for MethylResultSet
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges GRangesList
#' @importFrom dplyr left_join filter
setMethod("get_dmps_ranges", "MethylResultSet",
  function(x, index = 1, q = 0.05, keep_extra_columns = TRUE, as_list = FALSE) {
    if (is.null(index)) {
      warning(sprintf(
        "No index specified.",
        "Please provide a valid index or name.",
        "Returning empty GRanges."
      ))
      return(GenomicRanges::GRanges())
    }

    if (!index %in% names(x@dmps)) {
      warning(sprintf(
        "index '%s' not found in MethylResultSet",
        "contrasts. Returning empty GRanges.",
        index
      ))
      return(GenomicRanges::GRanges())
    }

    df <- dplyr::left_join(
      as.data.frame(x@manifest),
      x@dmps[[index]],
      by  = "Probe_ID"
    ) |> dplyr::filter(.data$adj.P.Val < q)

    if (!as_list) {
      return(GenomicRanges::makeGRangesFromDataFrame(
        df,
        seqnames.field = "CpG_chrm",
        start.field = "CpG_beg",
        end.field = "CpG_end",
        keep.extra.columns  = keep_extra_columns,
        na.rm = TRUE
      ))
    }

    # as_list = TRUE: split by deltabetas sign and all
    if (!"deltabetas" %in% colnames(df)) {
      stop("deltabetas column is required in DMP table for splitting.")
    }
    plateform <- x@metadata$plateform
    gr_all <- GenomicRanges::makeGRangesFromDataFrame(
      as.data.frame(x@manifest),
      seqnames.field = "CpG_chrm",
      start.field = "CpG_beg",
      end.field = "CpG_end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    gr_hyper <- GenomicRanges::makeGRangesFromDataFrame(
      dplyr::filter(df, .data$deltabetas > 0),
      seqnames.field = "CpG_chrm",
      start.field = "CpG_beg",
      end.field = "CpG_end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    gr_hypo <- GenomicRanges::makeGRangesFromDataFrame(
      dplyr::filter(df, .data$deltabetas < 0),
      seqnames.field = "CpG_chrm",
      start.field = "CpG_beg",
      end.field = "CpG_end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    grl <- GenomicRanges::GRangesList(
      hypermethylated = gr_hyper,
      hypomethylated = gr_hypo
    )
    grl[[plateform]] <- gr_all

    grl
  }
)



#' Retrieve DRP Ranges as a GRanges Object
#'
#' Returns significant DRPs for a given contrast as a \code{GRanges} object,
#' filtered by adjusted p-value.
#'
#' @param x A \code{MethylResultSet} object.
#' @param index Integer index or character name of the contrast
#' (default \code{1}).
#' @param q Adjusted p-value threshold (default \code{0.05}).
#'
#' @return A \code{GRanges} object with one range per significant probe.
#'
#' @export
setGeneric("get_dmrs_ranges",
  function(
    x,
    index,
    q  = 0.05,
    keep_extra_columns = TRUE,
    as_list = FALSE,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff")
  ) {
    standardGeneric("get_dmrs_ranges")
  }
)

#' @importFrom dplyr filter
#' @importFrom GenomicRanges makeGRangesFromDataFrame GRanges GRangesList
#' @describeIn get_dmrs_ranges Method for MethylResultSet
setMethod("get_dmrs_ranges", "MethylResultSet",
  function(
    x,
    index,
    q  = 0.05,
    keep_extra_columns = TRUE,
    as_list = FALSE,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff")
  ) {
    if (is.null(index)) {
      warning(sprintf(
        "No index specified.",
        "Please provide a valid index or name.",
        "Returning empty GRanges."
      ))
      return(GenomicRanges::GRanges())
    }

    if (!index %in% names(x@dmps)) {
      warning(sprintf(
        "index '%s' not found in MethylResultSet",
        "contrasts. Returning empty GRanges.",
        index
      ))
      return(GenomicRanges::GRanges())
    }

    df <- get_dmrs(x, index) |>
      dplyr::filter(.data$tool_fdr < q) |>
      dplyr::filter(.data$tool %in% tools)

    if (!as_list) {
      return(GenomicRanges::makeGRangesFromDataFrame(
        df,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns  = keep_extra_columns,
        na.rm = TRUE
      ))
    }

    # as_list = TRUE: split by deltabetas sign and all
    if (!"mean_deltabeta" %in% colnames(df)) {
      stop("mean_deltabeta column is required in DMP table for splitting.")
    }
    plateform <- x@metadata$plateform
    gr_all <- GenomicRanges::makeGRangesFromDataFrame(
      as.data.frame(x@manifest),
      seqnames.field = "CpG_chrm",
      start.field = "CpG_beg",
      end.field = "CpG_end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    gr_hyper <- GenomicRanges::makeGRangesFromDataFrame(
      dplyr::filter(df, .data$mean_deltabeta > 0),
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    gr_hypo <- GenomicRanges::makeGRangesFromDataFrame(
      dplyr::filter(df, .data$mean_deltabeta < 0),
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end",
      keep.extra.columns  = keep_extra_columns,
      na.rm = TRUE
    )
    grl <- GenomicRanges::GRangesList(
      hypermethylated = gr_hyper,
      hypomethylated = gr_hypo
    )
    grl[[plateform]] <- gr_all

    grl
  }
)

#' Retrieve Full Per-Probe Results with DMR Annotations
#'
#' Returns all probes for a given contrast annotated with their limma
#' statistics and DMR membership as reported by one or more DMR-calling tools.
#' The result is in \strong{long} format: one row per probe-DMR combination.
#'
#' @param x A \code{MethylResultSet} object.
#' @param index Integer index or character name of the contrast.
#' @param tools Character vector of DMR tools to run. Any subset of
#'   \code{c("dmrcate", "ipdmr", "combp", "dmrff")} (default
#'   \code{c("dmrcate", "ipdmr")}).
#' @param maxgap Maximum gap in base pairs between consecutive probes to be
#'   considered part of the same DMR (default \code{1000}).
#' @param mvals A numeric matrix of M-values (probes × samples). Required when
#'   \code{"dmrff"} is included in \code{tools}.
#' @param ... Additional arguments forwarded to \code{add_dmrcate}
#'   (e.g. \code{fdr}, \code{pcutoff}, \code{genome}).
#'
#' @return A data frame in long format with columns:
#'   \itemize{
#'     \item All columns from \code{\link{getDMPs}} output
#'     \item \code{dmrtool}: Tool name that detected the DMR
#'     \item \code{ID}: DMR identifier ("Chr-Start-End-tool")
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{fdr}: FDR value for the DMR
#'     \item \code{no.cpgs}: Number of CpGs in the DMR
#'     \item \code{tool}: Tool name (same as dmrtool)
#'   }
#'   Only probes that belong to at least one DMR are included.
#'
#' @export
setGeneric("get_results",
  function(
    x,
    index,
    tools = c("dmrcate", "ipdmr"),
    maxgap  = 1000,
    mvals  = NULL,
    genome  = "hg38", ...
  ) {
    standardGeneric("get_results")
  }
)

#' @importFrom dplyr bind_rows filter
#' @describeIn get_results Method for MethylResultSet
setMethod("get_results", "MethylResultSet",
  function(
    x,
    index,
    tools = c("dmrcate", "ipdmr"),
    maxgap = 1000,
    mvals = NULL,
    genome = "hg38", ...
  ) {

    if ("dmrff" %in% tools && is.null(mvals)) {
      stop("'mvals' must be provided when 'dmrff' is included in tools.")
    }
    dmps <- get_dmps(x, index) |>
      dplyr::filter(!is.na(.data$chr)) |>
      dplyr::filter(!is.na(.data$P.Value))

    # Collect DMR results from each tool
    dmr_results <- list()

    if ("dmrcate" %in% tools) {
      dmr_results <- c(dmr_results,
        list(add_dmrcate(dmps, genome = genome, ...))
      )
    }
    if ("ipdmr" %in% tools) {
      dmr_results <- c(dmr_results, list(add_ipdmr(dmps, maxgap)))
    }
    if ("combp" %in% tools) {
      dmr_results <- c(dmr_results, list(add_combp(dmps, maxgap)))
    }
    if ("dmrff" %in% tools) {
      dmr_results <- c(dmr_results, list(add_dmrff(dmps, mvals, maxgap)))
    }

    if (length(dmr_results) == 0) return(data.frame())

    dmr_combined <- dplyr::bind_rows(dmr_results)

    x@dmrs[[index]] <- dmr_combined
    x
  }
)

#' Add dmrcate DMR Annotations to a DMP Table
#'
#' Runs \code{dmrtools_dmrcate} on the supplied DMP data frame and returns
#' DMR annotations in long format with one row per probe-DMR combination.
#'
#' @param x A data frame of DMPs as returned by \code{\link{getDMPs}} (must
#'   contain at least \code{Probe_ID} and \code{P.Value}).
#' @param fdr FDR threshold for DMR calling (default \code{0.05}).
#' @param pcutoff Probe-level p-value cut-off passed to dmrcate (default
#'   \code{0.05}).
#' @param maxgap Maximum gap between probes in a DMR, in bp (default
#'   \code{1000}).
#' @param genome Genome assembly identifier (default \code{"hg38"}).
#' @param mean Statistic used to summarise DMR-level significance. One of
#'   \code{"HMFDR"} (default) or \code{"min_smoothed_fdr"}.
#'
#' @return A data frame in long format with columns:
#'   \itemize{
#'     \item \code{dmrtool}: Tool name ("dmrcate")
#'     \item \code{ID}: DMR identifier ("Chr-Start-End-dmrcate")
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{fdr}: FDR value
#'     \item \code{no.cpgs}: Number of CpGs in DMR
#'     \item \code{tool}: Tool name (same as dmrtool)
#'     \item \code{Probe_ID}: Probe identifier
#'   }
#'   Returns an empty data frame if no DMRs are found.
#'
#' @export
setGeneric("add_dmrcate",
  function(
    x,
    fdr  = 0.2,
    pcutoff  = 0.05,
    maxgap  = 1000,
    genome  = "hg38",
    mean  = "HMFDR"
  ) {
    standardGeneric("add_dmrcate")
  }
)

#' @importFrom dplyr mutate select
#' @importFrom tidyr separate_longer_delim
#' @describeIn add_dmrcate Method for MethylResultSet
setMethod("add_dmrcate", "data.frame",
  function(
    x,
    fdr  = 0.2,
    pcutoff  = 0.05,
    maxgap  = 1000,
    genome  = "hg38",
    mean  = "HMFDR"
  ) {

    if (!mean %in% c("HMFDR", "min_smoothed_fdr")) {
      warning(
        "parameter mean should be 'HMFDR' or 'min_smoothed_fdr';",
        " defaulting to HMFDR"
      )
      mean <- "HMFDR"
    }

    dmrs <- dmrtools_dmrcate(
      x,
      fdr  = fdr,
      pcutoff = pcutoff,
      maxgap = maxgap,
      genome  = genome
    )

    if (is.null(dmrs) || nrow(dmrs) == 0) {
      return(data.frame())
    }

    dmrs <- dmrs |>
      tidyr::separate_longer_delim("Probe_ID", delim = ";") |>
      dplyr::mutate(
        dmrtool = "dmrcate",
        ID = paste(seqnames, start, end, "dmrcate", sep = "-"),
        tool = "dmrcate"
      ) |>
      dplyr::select(
        "dmrtool",
        "ID",
        Start = "start",
        End = "end",
        fdr = !!mean,
        "no.cpgs",
        "tool",
        "Probe_ID"
      )

    dmrs
  }
)

#' Add ipdmr DMR Annotations to a DMP Table
#'
#' Runs \code{dmrtools_ipdmr} on the supplied DMP data frame and returns
#' DMR annotations in long format with one row per probe-DMR combination.
#'
#' @param x A data frame of DMPs (must contain \code{Probe_ID}, genomic
#'   coordinates, and \code{P.Value}).
#' @param maxgap Maximum gap between probes in a DMR, in bp (default
#'   \code{1000}).
#'
#' @return A data frame in long format with columns:
#'   \itemize{
#'     \item \code{dmrtool}: Tool name ("ipdmr")
#'     \item \code{ID}: DMR identifier ("Chr-Start-End-ipdmr")
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{fdr}: FDR value
#'     \item \code{no.cpgs}: Number of CpGs in DMR
#'     \item \code{tool}: Tool name (same as dmrtool)
#'     \item \code{Probe_ID}: Probe identifier
#'   }
#'   Returns an empty data frame if no DMRs are found.
#'
#' @export
setGeneric("add_ipdmr", function(x, maxgap  = 1000) {
  standardGeneric("add_ipdmr")
})

#' @importFrom dplyr mutate select
#' @importFrom tidyr separate_longer_delim
#' @importFrom stringr str_starts
#' @describeIn add_ipdmr Method for MethylResultSet
setMethod(
  "add_ipdmr", "data.frame",
  function(x, maxgap  = 1000) {
    dmrs <- dmrtools_ipdmr(x, maxgap = maxgap)
    if (is.null(dmrs) || nrow(dmrs) == 0) {
      return(data.frame())
    }

    dmrs <- dmrs |>
      tidyr::separate_longer_delim("probe", delim = ";") |>
      dplyr::rename(Probe_ID = "probe") |>
      dplyr::mutate(chr = as.character(.data$chr)) |>
      dplyr::mutate(
        dmrtool = "ipdmr",
        chr = ifelse(
          !stringr::str_starts(.data$chr, "chr"),
          paste0("chr", .data$chr),
          .data$chr
        ),
        ID = paste(.data$chr, .data$start, .data$end, "ipdmr", sep = "-"),
        tool = "ipdmr"
      ) |>
      dplyr::select(
        "dmrtool",
        "ID",
        Start = "start",
        End = "end",
        "fdr",
        no.cpgs = "nprobe",
        "tool",
        "Probe_ID"
      )

    dmrs

  }
)

#' Add comb-p DMR Annotations to a DMP Table
#'
#' Runs \code{dmrtools_combp} on the supplied DMP data frame and returns
#' DMR annotations in long format with one row per probe-DMR combination.
#'
#' @param x A data frame of DMPs (must contain \code{Probe_ID}, genomic
#'   coordinates, and \code{P.Value}).
#' @param maxgap Maximum gap between probes in a DMR, in bp (default
#'   \code{1000}).
#'
#' @return A data frame in long format with columns:
#'   \itemize{
#'     \item \code{dmrtool}: Tool name ("combp")
#'     \item \code{ID}: DMR identifier ("Chr-Start-End-combp")
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{fdr}: FDR value
#'     \item \code{no.cpgs}: Number of CpGs in DMR
#'     \item \code{tool}: Tool name (same as dmrtool)
#'     \item \code{Probe_ID}: Probe identifier
#'   }
#'   Returns an empty data frame if no DMRs are found.
#'
#' @export
setGeneric("add_combp", function(x, maxgap  = 1000) {
  standardGeneric("add_combp")
})

#' @importFrom dplyr mutate select
#' @importFrom tidyr separate_longer_delim
#' @importFrom stringr str_starts
#' @describeIn add_combp Method for MethylResultSet
setMethod("add_combp", "data.frame", function(x, maxgap  = 1000) {
  dmrs <- dmrtools_combp(x, maxgap  = maxgap)
  if (is.null(dmrs) || nrow(dmrs) == 0) {
    return(data.frame())
  }

  dmrs <- dmrs |>
    tidyr::separate_longer_delim("probe", delim = ";") |>
    dplyr::rename(Probe_ID = "probe") |>
    dplyr::mutate(
      dmrtool = "combp",
      chr = ifelse(
        !stringr::str_starts(.data$chr, "chr"),
        paste0("chr", .data$chr),
        .data$chr
      ),
      ID = paste(.data$chr, .data$start, .data$end, "combp", sep = "-"),
      tool = "combp"
    ) |>
    dplyr::select(
      "dmrtool",
      "ID",
      Start = "start",
      End = "end",
      "fdr",
      no.cpgs = "nprobe",
      "tool",
      "Probe_ID"
    )

  dmrs
})

#' Add dmrff DMR Annotations to a DMP Table
#'
#' Runs \code{dmrtools_dmrff} on the supplied DMP data frame and appends a
#' \code{dmrff} column encoding DMR membership as
#' \code{"chr:start-end:fdr:nprobe"}.
#'
#' @param x A data frame of DMPs (must contain \code{Probe_ID}, genomic
#'   coordinates, and \code{P.Value}).
#' @param mvals A numeric matrix of M-values (probes × samples). Required by
#'   the dmrff algorithm.
#' @param maxgap Maximum gap between probes in a DMR, in bp (default
#'   \code{1000}).
#'
#' @return The input data frame with an additional \code{dmrff} column.
#'   Returns \code{x} unchanged if no DMRs are found.
#'
#' @export
setGeneric("add_dmrff", function(x, mvals, maxgap  = 1000) {
  standardGeneric("add_dmrff")
})

#' @importFrom dplyr mutate select
#' @importFrom tidyr separate_longer_delim
#' @importFrom stringr str_starts
#' @describeIn add_dmrff Method for MethylResultSet
setMethod("add_dmrff", "data.frame", function(x, mvals, maxgap  = 1000) {
  dmrs <- dmrtools_dmrff(x, mvals, maxgap  = maxgap)
  if (is.null(dmrs) || nrow(dmrs) == 0) {
    return(data.frame())
  }

  dmrs <- dmrs |>
    tidyr::separate_longer_delim("Probe_ID", delim = ";") |>
    dplyr::mutate(
      dmrtool = "dmrff",
      chr = ifelse(
        !stringr::str_starts(.data$chr, "chr"),
        paste0("chr", .data$chr),
        .data$chr
      ),
      ID = paste(.data$chr, .data$start, .data$end, "dmrff", sep = "-"),
      tool = "dmrff"
    ) |>
    dplyr::select(
      "dmrtool",
      "ID",
      Start = "start",
      End = "end",
      "fdr",
      no.cpgs = "nprobe",
      "tool",
      "Probe_ID"
    )

  dmrs
})


#' Retrieve Full Per-Probe Results with DMR Annotations
#'
#' Returns all probes for a given contrast annotated with their limma
#' statistics and DMR membership as reported by one or more DMR-calling tools.
#' The result is in \strong{long} format: one row per probe-DMR combination.
#'
#' @param mrs A \code{MethylResultSet} object.
#' @param index Integer index or character name of the contrast.
#'
#' @return A data frame in long format with columns:
#'   \itemize{
#'     \item All columns from \code{\link{getDMPs}} output
#'     \item \code{dmrtool}: Tool name that detected the DMR
#'     \item \code{ID}: DMR identifier ("Chr-Start-End-tool format")
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{fdr}: FDR value for the DMR
#'     \item \code{no.cpgs}: Number of CpGs in the  DMR
#'     \item \code{tool}: Tool name (same as dmrtool)
#'     \item \code{Probe_ID}: Probe identifier
#'     \item \code{chr}: Chromosome (from manifest)
#'     \item \code{pos}: Genomic position (from manifest)
#'   }
#'   All dmps and combinations of dmps and dmrs are included.
#'
#' @importFrom dplyr left_join arrange
#'
#' @export
setGeneric("as_dataframe",
  function(
    mrs,
    index
  ) {
    standardGeneric("as_dataframe")
  }
)

#' @importFrom dplyr left_join arrange
#'
#' @describeIn as_dataframe Method for MethylResultSet
setMethod("as_dataframe", "MethylResultSet",
  function(
    mrs,
    index
  ) {

    dplyr::left_join(
      as.data.frame(mrs@dmps[[index]]),
      as.data.frame(mrs@dmrs[[index]]),
      by = "Probe_ID"
    ) |>
      dplyr::left_join(
        as.data.frame(mrs@manifest),
        by = "Probe_ID"
      ) |>
      dplyr::arrange(.data$adj.P.Val, .data$fdr, .data$no.cpgs)
  }
)

#' Get Unique DMR Table with Summary Statistics
#'
#' Aggregates DMR results from \code{\link{get_results}}
#' into a table of unique DMRs with summary statistics and probe lists.
#'
#' @param results A data frame as returned by
#'   \code{\link{get_results}} (long format).
#' @param tools Character vector of DMR tools to include.
#'   If \code{NULL} (default), includes all tools present in the results.
#' @param max_fdr Maximum FDR threshold for DMR inclusion (default \code{0.05}).
#' @param min_cpgs Minimum number of CpGs required for DMR inclusion
#'   (default \code{2}).
#'
#' @return A data frame with one row per unique DMR, containing:
#'   \itemize{
#'     \item \code{ID}: Unique DMR identifier
#'     \item \code{chr}: Chromosome
#'     \item \code{Start}: DMR start position
#'     \item \code{End}: DMR end position
#'     \item \code{tools}: Comma-separated list of tools that detected this DMR
#'     \item \code{mean_fdr}: Mean FDR across detecting tools
#'     \item \code{min_fdr}: Minimum FDR across detecting tools
#'     \item \code{no.cpgs}: Number of CpGs in the DMR
#'     \item \code{probes}: Semicolon-separated list of probe IDs in the DMR
#'     \item \code{mean_deltabeta}: Arithmetic mean of delta beta values
#'     \item \code{geometric_mean_deltabeta}:
#'       Geometric mean of delta beta values
#'     \item \code{harmonic_mean_deltabeta}:
#'       Harmonic mean of delta beta values
#'   }
#'
#' @details
#' Delta beta averages are calculated using three different methods:
#' \itemize{
#'   \item Arithmetic mean: \code{mean(deltabetas)}
#'   \item Geometric mean: \code{exp(mean(log(abs(deltabetas) + 1e-10))) *
#'     sign(mean(deltabetas))}
#'   \item Harmonic mean: \code{1/mean(1/(abs(deltabetas) + 1e-10)) * #'     sign(mean(deltabetas))}
#' }
#'
#' @export
setGeneric("get_dmrs", function(
    mrs,
    index,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
    max_fdr = 0.05,
    min_cpgs = 2) {
  standardGeneric("get_dmrs")
})

#' @importFrom dplyr group_by summarize filter arrange
#' @describeIn get_dmrs Method for MethylResultSet
setMethod("get_dmrs", "MethylResultSet", function(
    mrs,
    index,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
    max_fdr = 0.05,
    min_cpgs = 2) {

  list_uniq <- function(column) {
    column |>
      paste(collapse = ";") |>
      strsplit(";") |>
      unlist() |>
      trimws() |>
      (\(x) x[x != "NA"])() |>
      unique() |>
      paste(collapse = ";")
  }

  # Group by DMR ID and summarize
  as_dataframe(mrs, index) |>
    dplyr::filter(.data$tool %in% tools) |>
    dplyr::filter(.data$fdr < max_fdr) |>
    dplyr::filter(.data$no.cpgs >= min_cpgs) |>
    dplyr::group_by(.data$ID, .data$tool) |>
    dplyr::rename(any_of(c(Feature_UCSC = "UCSC_RefGene_Group"))) |>
    dplyr::summarize(
      chr = dplyr::first(.data$CpG_chrm),
      start = min(.data$Start),
      end = max(.data$End),
      tools = paste(sort(unique(.data$tool)), collapse = ","),
      tool_fdr = dplyr::first(.data$fdr),
      min_fdr = min(.data$adj.P.Val, na.rm = TRUE),
      max_fdr = max(.data$adj.P.Val, na.rm = TRUE),
      HMFDR = h_mean(.data$adj.P.Val, na.rm = TRUE),
      no.cpgs = dplyr::first(.data$no.cpgs),
      probes = paste(sort(unique(.data$Probe_ID)), collapse = ";"),
      mean_deltabeta = mean(.data$deltabetas, na.rm = TRUE),
      mean_abs_deltabeta = mean(abs(.data$deltabetas), na.rm = TRUE),
      max_deltabeta = max(.data$deltabetas, na.rm = TRUE),
      Relation_to_Island = list_uniq(.data$Relation_to_Island),
      Feature_UCSC = list_uniq(.data$Feature_UCSC),
      genesUniq = list_uniq(.data$genesUniq),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$chr, .data$start)

})

#' Extract Top Probes with Extended Statistics
#'
#' Extracts the top probes for a given contrast from a \code{limma}
#'   eBayes object,
#' and appends additional per-probe statistics including regression
#'  coefficients,
#' standard deviations, and goodness-of-fit (R-squared) values.
#'
#' @param eb An eBayes object as returned by \code{limma::eBayes}.
#' @param x Coefficient index or name specifying which contrast to extract
#'   (passed to \code{limma::topTable}).
#' @param rsq A numeric vector of goodness-of-fit (R-squared) values,
#'   one per probe.
#'   Must have names matching probe IDs.
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item All columns from \code{limma::topTable} output
#'     \item \code{Coefficient}: Regression coefficient for the
#'       selected contrast
#'     \item \code{Stdev}: Posterior standard deviation (shrunk estimate)
#'     \item \code{goodness}: R-squared (goodness of fit) for each probe
#'   }
#'
#' @details
#' This is an internal utility function used by the \code{MethylResultSet}
#' constructor to annotate DMPs with additional model fit statistics.
#'
#' @importFrom limma topTable
top_tables <- function(eb, x, rsq) {

  x_table <- limma::topTable(
    eb,
    adjust = "BH",
    number = Inf,
    p = 1,
    sort.by = "P",
    coef = x
  )

  x_table$Coefficient <- eb$coefficients[rownames(x_table), x]
  x_table$Stdev <- (sqrt(eb$s2.post) * eb$stdev.unscaled)[rownames(x_table), x]
  x_table$goodness <- rsq[rownames(x_table)]

  x_table
}


#' Generate Model Name from MethylResultSet Metadata
#'
#' Constructs a descriptive, human-readable name for a model analysis
#' by combining metadata components (SVA correction,
#'   model formula, fitting method,
#' and intercept level) in a standardized format suitable for file naming.
#'
#' @param x A \code{MethylResultSet} object.
#' @param format Format style for the model name. One of:
#'   \itemize{
#'     \item \code{"clean"} (default): Clean format with descriptive labels.
#'       Example: \code{sva-group_ls_WT_mod-group-age-sex}
#'     \item \code{"abbreviated"}: Abbreviated format with shortened components.
#'       Example: \code{s1_ls_WT_mod-grp-age-gen-id}
#'   }
#'
#' @return A character string containing the formatted model name.
#'
#' @details
#' The function extracts from \code{x@metadata}:
#' \itemize{
#'   \item \code{sva}: SVA correction status and protected variables
#'   \item \code{model}: Model formula used for fitting
#'   \item \code{method}: Fitting method ("ls", "robust", etc.)
#'   \item \code{intercept}: Reference level for the grouping factor
#' }
#'
#' \strong{Clean format}:
#' \itemize{
#'   \item SVA: \code{sva-0} (no correction) or \code{sva-adj-age-sex} 
#'     (protected variables)
#'   \item Method: \code{ls} or \code{rl}
#'   \item Intercept: Reference level (e.g., \code{WT})
#'   \item Model: Variables separated by hyphens (mathematical symbols replaced)
#' }
#'
#' \strong{Abbreviated format}:
#' \itemize{
#'   \item SVA: \code{S0} (no SVA) or \code{S1} (SVA applied)
#'   \item Method: \code{ls} or \code{rl}
#'   \item Intercept: Reference level (e.g., \code{WT})
#'   \item Model: First 3 letters of each variable
#' }
#'
#' @examples
#' # Assuming \code{mrs} is a MethylResultSet object
#' \dontrun{
#'   get_model_name(mrs)  # Clean format (default)
#'   get_model_name(mrs, format = "abbreviated")
#' }
#'
#' @export
setGeneric("get_model_name", function(x, format = "clean"){
  standardGeneric("get_model_name")
})

#' @describeIn get_model_name Method for MethylResultSet
setMethod("get_model_name", "MethylResultSet", function(x, format = "clean") {

  # Extract metadata components
  sva        <- x@metadata$sva %||% "unknown"
  model      <- x@metadata$model %||% "unknown"
  method     <- x@metadata$method %||% "ls"
  intercept  <- x@metadata$intercept %||% "ref"

  # Format based on user choice
  if (format == "clean") {
    # Clean format: sva-group_ls_WT_mod-group-age-sex

    # Process SVA component
    if (is.null(sva) || sva == "no" || sva == "unknown") {
      sva_part <- "sva-0"
    } else if (sva == "yes") {
      sva_part <- "sva-adj"
    } else {
      # Assume sva contains the protected variables list
      sva_part <- paste0("sva-", sva)
    }

    # Process model component: replace mathematical symbols with hyphens
    model_clean <- gsub("[~+\\*\\(\\)]", "-", model)
    # Remove multiple consecutive hyphens
    model_clean <- gsub("-+", "-", model_clean)
    # Remove any underscores in variable names to avoid confusion
    model_clean <- gsub("_", "", model_clean)
    # Remove leading/trailing hyphens
    model_clean <- trimws(model_clean, whitespace = "-")

    # Method component (standardize to ls/rl)
    method_part <- ifelse(method == "robust", "rl", "ls")

    # Combine components: [SVA][Méthode]_[Intercept]_[Modèle]
    name <- sprintf(
      "%s_%s_%s_mod-%s",
      sva_part, method_part, intercept, model_clean
    )

  } else if (format == "abbreviated") {
    # Abbreviated format: s1_ls_WT_mod-grp-age-gen-id

    # Process SVA component
    if (is.null(sva) || sva == "no" || sva == "unknown") {
      sva_part <- "s0"
    } else {
      sva_part <- "s1"
    }

    # Process model component: take first 3 letters of each variable
    # Split by mathematical operators and take first 3 chars of each term
    terms <- unlist(strsplit(model, "[~+\\*\\(\\)]"))
    terms <- terms[terms != ""]  # Remove empty strings
    abbr_terms <- sapply(terms, function(term) {
      substr(trimws(term), 1, 3)
    })
    model_abbr <- paste(abbr_terms, collapse = "-")

    # Method component (standardize to ls/rl)
    method_part <- ifelse(method == "robust", "rl", "ls")

    # Combine components: [SVA][Méthode]_[Intercept]_[Modèle]
    name <- sprintf(
      "%s_%s_%s_mod-%s",
      sva_part, method_part, intercept, model_abbr
    )

  } else {
    stop("format must be 'clean' or 'abbreviated'")
  }

  name
})

#' Count Unique DMRs Detected by Each Tool
#'
#' Aggregates DMR results to count the number of unique DMRs detected by each
#' DMR-calling tool.
#'
#' @param df A data frame containing DMR results in long format, as returned by
#'   \code{\link{get_results}}. Must contain columns \code{ID}
#'   (unique DMR identifier) and \code{dmrtool}
#'   (name of the DMR tool that detected the DMR).
#'   The function will count unique DMRs based on the \code{ID}
#'   and \code{dmrtool}.
#'
#' @return A named integer vector where names are DMR tools and values are the
#'   count of unique DMRs detected by each tool.
#' @export
setGeneric("get_number_of_dmrs", function(x, index) {
  standardGeneric("get_number_of_dmrs")
})


#' @describeIn get_number_of_dmrs Method for MethylResultSet
setMethod("get_number_of_dmrs", "MethylResultSet",
  function(x, index) {

    if (length(x@dmrs[[index]]) == 0) {
      return(setNames(integer(0), character(0)))
    }

    x@dmrs[[index]] |>
      select("ID", "tool") |>
      unique() |>
      pull("tool") |>
      table()
  }
)


#' Extract Unique Genes from Significant DMPs and DMRs
#' Retrieves a list of unique gene symbols associated with significant DMPs and
#' DMRs for a given contrast, based on specified significance thresholds and DMR
#' tools.
#' @param x A \code{MethylResultSet} object containing DMP and DMR results.
#' @param index Integer index or character name of the contrast to analyze.
#' @param tools Character vector of DMR tools to consider for DMR-based gene
#'   extraction. Default is \code{c("dmrcate", "ipdmr", "combp", "dmrff")}.
#' @param q Numeric value specifying the FDR threshold for significance
#'   (default \code{0.05}).
#' @param deltabetas Numeric value specifying the minimum absolute
#'   delta beta threshold for significance (default \code{0.3}).
#' @param min_cpgs Integer specifying the minimum number of CpGs
#'   required in a DMR for it to be considered significant (default \code{2}).
#' @return A character vector of unique gene symbols associated with
#'   significant DMPs and DMRs based on the specified criteria.
#' @details The function performs the following steps:
#' \itemize{
#'   \item Extracts significant DMPs based on the specified FDR and
#'     delta beta thresholds, and collects associated gene symbols.
#'   \item Extracts significant DMRs based on the specified DMR tools,
#'     FDR threshold, delta beta threshold, and minimum CpG count,
#'     and collects associated gene symbols.
#'   \item Combines gene symbols from both DMPs and DMRs, removes duplicates,
#'     and returns a unique list of gene symbols.
#' }
#' @export
setGeneric("get_genes",
  function(
    x,
    index,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
    q = 0.05,
    deltabetas = 3.0,
    min_cpgs = 2
  ) {
    standardGeneric("get_genes")
  }
)

#' @importFrom dplyr filter select distinct pull
#' @describeIn get_genes Method for MethylResultSet
setMethod("get_genes", "MethylResultSet",
  function(
    x,
    index,
    tools = c("dmrcate", "ipdmr", "combp", "dmrff"),
    q = 0.05,
    deltabetas = 3.0,
    min_cpgs = 2
  ) {

    # genes from significant dmps
    dmps_genes <- get_dmps(x, index) |>
      dplyr::filter(!is.na(.data$genesUniq)) |>
      dplyr::filter(.data$adj.P.Val < q) |>
      dplyr::filter(.data$deltabetas >= deltabetas)

    # genes from significant dmrs
    dmrs_genes <- get_dmrs(x, index) |>
      dplyr::filter(!is.na(.data$genesUniq)) |>
      dplyr::filter(.data$tool %in% tools) |>
      dplyr::filter(.data$tool_fdr < q) |>
      dplyr::filter(.data$mean_abs_deltabeta >= deltabetas) |>
      dplyr::filter(.data$no.cpgs >= min_cpgs)

    dplyr::bind_rows(dmps_genes, dmrs_genes) |>
      dplyr::select("genesUniq") |>
      dplyr::distinct() |>
      dplyr::pull("genesUniq") |>
      strsplit(";") |>
      unlist() |>
      trimws() |>
      unique()

  }
)
