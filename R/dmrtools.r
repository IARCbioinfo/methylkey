#' Search for Differentially Methylated Regions (DMRs) using DMRcate
#'
#' This function searches for Differentially Methylated Regions (DMRs)
#'   in DNA methylation data using the DMRcate package.
#'   It takes a set of CpG sites with associated statistical information
#'   and annotates them for DMR analysis.
#'
#' @param dmps A data frame containing CpG site information,
#'   including chromosome, position, strand, statistical test statistic,
#'   delta-betas, adjusted p-values, and probe IDs.
#' @param fdr The q-value threshold for filtering CpG sites to be considered
#'   significant. Default is 0.2.
#' @param maxgap The maximum gap between neighboring CpGs to consider
#'   when defining DMRs. Default is 1000.
#' @param pcutoff Threshold to determine DMRs. Default is 0.05.
#' @param genome The genome build to use for genomic annotation.
#'   Default is "hg19."
#'
#' @return A data frame with information about identified DMRs,
#'   including chromosome, start and end positions,
#'   and associated CpG probe IDs.
#'
#' @importFrom dplyr filter bind_cols
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom DMRcate dmrcate extractRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stringr str_detect
#'
#' @export
dmrtools_dmrcate <- function(
    dmps,
    fdr = 0.2,
    maxgap = 1000,
    pcutoff = 0.05,
    genome = "hg38") {

  if (!requireNamespace("DMRcate", quietly = TRUE)) {
    stop("Package 'DMRcate' is required for this function to work",
         "Please install it.",
         call. = FALSE)
  }

  dmps <- dmps |> dplyr::filter(!stringr::str_detect(.data$chr, "_"))

  annotated <- data.frame(
    chr = dmps$chr,
    start = dmps$pos,
    end = dmps$pos,
    strand = dmps$strand,
    rawpval = dmps$P.Value,
    stat = dmps$t,
    diff = dmps$deltabetas,
    ind.fdr = dmps$adj.P.Val,
    is.sig = (dmps$adj.P.Val < fdr)
  )
  annotated <- GenomicRanges::makeGRangesFromDataFrame(
    annotated,
    keep.extra.columns = TRUE
  )
  names(annotated) <- dmps$Probe_ID
  myannotation <- new("CpGannotated", ranges = sort(annotated))
  if (sum(is.na(myannotation@ranges$diff))) {
    myannotation@ranges$diff[which(is.na(myannotation@ranges$diff))] <- 0
  }
  if (sum(annotated$is.sig) < 1) {
    warning("No valid CpG sites remaining after filtering")
    return(NULL)
  }

  dmrcoutput <- DMRcate::dmrcate(
    myannotation,
    C = 2,
    pcutoff = pcutoff,
    lambda = maxgap
  )
  table <- DMRcate::extractRanges(dmrcoutput, genome = genome)

  overlap <- GenomicRanges::findOverlaps(annotated, table, type = "within")
  table <- as.data.frame(table)[S4Vectors::subjectHits(overlap), c(
    "seqnames", "start", "end",
    "HMFDR", "min_smoothed_fdr", "no.cpgs"
  )]

  dmps[S4Vectors::queryHits(overlap), ] |> dplyr::bind_cols(table)
}

#' Search for Differentially Methylated Regions (DMRs) using DMRff
#'
#' This function searches for Differentially Methylated Regions (DMRs)
#'   in DNA methylation data using the DMRff package. It takes a set of
#'   CpG sites with associated statistical information and methylation data
#'   and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information, including
#'   chromosome, position, statistical coefficient, standard deviation,
#'   p-value, and probe IDs.
#' @param betas A data frame containing methylation data, typically beta values,
#'   for the CpG sites.
#' @param maxgap The maximum gap between neighboring CpGs to consider when
#'   defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs,
#'   including chromosome, start and end positions, associated CpG probe IDs,
#'   and DMR tool name.
#'
#'  @importFrom dplyr select mutate bind_cols
#'  @importFrom ENmix dmrff dmrff.sites
#'  @importFrom stats p.adjust

#' @export
dmrtools_dmrff <- function(dmps, betas, maxgap = 1000) {

  if (!requireNamespace("dmrff", quietly = TRUE)) {
    stop("Package 'dmrff' is required for this function to work.",
         "Please install it.",
         call. = FALSE)
  }

  dmrs <- dmrff::dmrff(
    estimate = as.vector(dmps$Coefficient),
    se = as.vector(dmps$Stdev),
    p.value = as.vector(dmps$P.Value),
    methylation = betas[dmps$Probe_ID, ],
    chr = as.vector(dmps$chr),
    pos = as.vector(dmps$pos),
    maxgap = maxgap,
    verbose = TRUE
  )

  dmrs <- as.data.frame(dmrs)
  dmrs <- dmrs[(dmrs$n >= 2), ]
  if (nrow(dmrs) < 1) {
    return(NULL)
  }
  dmrs$fdr <- stats::p.adjust(dmrs$p.value, method = "fdr")
  dmrs <- dmrs |>
    dplyr::mutate(
      ID = paste0(.data$chr, ":", .data$start, "-", .data$end),
      nprobe = .data$n,
      p = .data$p.value
    )

  sites <- dmrff::dmrff.sites(dmrs, dmps$chr, dmps$pos)
  final <- dplyr::bind_cols(
    Probe_ID = dmps[sites$site, "Probe_ID"],
    dmrs[sites$region, ]
  ) |>
    dplyr::select(.data$Probe_ID,
                  .data$chr,
                  .data$start,
                  .data$end,
                  .data$fdr,
                  .data$nprobe,
                  .data$p)

  final
}

#' Search for Differentially Methylated Regions (DMRs) using combp
#'
#' This function searches for Differentially Methylated Regions (DMRs)
#'   in DNA methylation data using the combp package.
#'   It takes a set of CpG sites with associated statistical information
#'   and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information,
#'   including chromosome, position, p-value, and probe IDs.
#' @param maxgap The maximum gap between neighboring CpGs to consider
#'   when defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs,
#'   including chromosome, start and end positions,
#'   associated CpG probe IDs, and DMR tool name.
#'
#' @importFrom readr write_csv read_csv
#' @importFrom dplyr filter select mutate
#' @importFrom ENmix combp
#' @importFrom stats p.adjust

#' @export
dmrtools_combp <- function(dmps, maxgap = 1000) {

  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required for this function to work",
         "Please install it.",
         call. = FALSE)
  }

  data <- data.frame(
    probe = dmps$Probe_ID,
    p = dmps$P.Value,
    chr = dmps$chr,
    start = dmps$pos,
    end = dmps$pos
  )

  # Generate empty table results to avoid to reload
  # previous results if combp find 0 dmrs.
  combp <- data.frame(
    chr = character(),
    start = numeric(),
    end = numeric(),
    p = numeric(),
    fdr = numeric(),
    nprobe = numeric(),
    probe = character()
  )
  readr::write_csv(combp, "resu_combp.csv")

  ENmix::combp(data,
    dist.cutoff = 1000,
    bin.size = 310,
    seed = 0.05,
    region_plot = FALSE,
    mht_plot = FALSE,
    nCores = 1,
    verbose = TRUE
  )
  dmrs <- readr::read_csv("resu_combp.csv")
  dmrs <- dmrs |>
    dplyr::filter(.data$nprobe >= 2) |>
    dplyr::select(-dplyr::all_of("fdr"))
  dmrs <- dmrs |>
    dplyr::mutate(fdr = stats::p.adjust(dmrs$p, method = "fdr"))
}

#' Search for Differentially Methylated Regions (DMRs) using ipdmr
#'
#' This function searches for Differentially Methylated Regions (DMRs)
#'   in DNA methylation data using the ipdmr package.
#'   It takes a set of CpG sites with associated statistical
#'   information and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information,
#'   including chromosome, position, p-value, and probe IDs.
#' @param maxgap The maximum gap between neighboring CpGs to consider
#'   when defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs,
#'   including chromosome, start and end positions, associated CpG
#'   probe IDs, and DMR tool name.
#'
#' @importFrom readr write_csv read_csv
#' @importFrom dplyr filter select mutate
#' @importFrom ENmix ipdmr
#' @importFrom stats p.adjust
#'
#' @export
dmrtools_ipdmr <- function(dmps, maxgap = 1000, bin_size = 310, seed = 0.05) {

  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required for this function to work",
         "Please install it.",
         call. = FALSE)
  }

  data <- data.frame(
    probe = dmps$Probe_ID,
    p = dmps$P.Value,
    chr = dmps$chr,
    start = dmps$pos,
    end = dmps$pos
  )
  data$p[is.na(data$p)] <- 0.99999

  # Generate empty table results to avoid to reload
  # previous results if combp find 0 dmrs.
  ipdmr <- data.frame(
    chr = character(),
    start = numeric(),
    end = numeric(),
    p = numeric(),
    fdr = numeric(),
    nprobe = numeric(),
    probe = character()
  )
  write_csv(ipdmr, "resu_ipdmr.csv")

  ENmix::ipdmr(data,
    dist.cutoff = maxgap,
    bin.size = bin_size,
    seed = seed,
    region_plot = FALSE,
    mht_plot = FALSE,
    verbose = FALSE
  )
  dmrs <- readr::read_csv("resu_ipdmr.csv") |>
    dplyr::filter(.data$nprobe > 1)

  dmrs
}
