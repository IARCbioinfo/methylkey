#' circosplot
#'
#' Publication-ready circos plot of DMPs/DMRs.
#'
#' @param ranges GRanges object containing at least `deltabetas`
#'   and optionally `genesUniq`.
#' @param genome Genome build: `hg38`, `hg19` or `mm10`.
#' @param label_probes Optional vector of probe IDs to annotate.
#'   If NULL and max_labels > 0, labels are chosen automatically.
#' @param max_labels Maximum number of labels to display
#'   when `label_probes` is NULL.
#' @param label_genes Optional character vector of gene symbols
#'   to annotate on the plot.
#' @param label_column Column containing gene annotations.
#'   Default is `genesUniq`.
#' @param label_size Text size for gene labels.
#' @param point_size Point size for circos dots.
#'
#' @return ggplot object
#'
#' @importFrom ggbio ggbio circle
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom Seqinfo seqlengths seqlevels seqinfo
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @export
#'
circosplot <- function(ranges, genome, label_probes = NULL,
                       label_size = 2.5, point_size = 1.2) {

  if (genome == "hg38") {
    require(BSgenome.Hsapiens.UCSC.hg38)
    species <- BSgenome.Hsapiens.UCSC.hg38
  } else if (genome == "hg19") {
    require(BSgenome.Hsapiens.UCSC.hg19)
    species <- BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "mm10") {
    require(BSgenome.Mmusculus.UCSC.mm10)
    species <- BSgenome.Mmusculus.UCSC.mm10
  } else {
    stop("Unsupported genome. Use hg38, hg19, or mm10.")
  }

  if (!"deltabetas" %in% names(mcols(ranges))) {
    stop("The GRanges object must contain a 'deltabetas' metadata column.")
  }

  if (!"midpoint" %in% names(mcols(ranges))) {
    ranges$midpoint <- start(ranges) + width(ranges) %/% 2
  }

  if (!"Probe_ID" %in% names(mcols(ranges))) {
    ranges$Probe_ID <- names(ranges)
    if (is.null(names(ranges)) || all(is.na(names(ranges)))) {
      ranges$Probe_ID <- paste0("probe_", seq_along(ranges))
    }
  }

  chr_len <- Seqinfo::seqlengths(species)
  chr_len <- chr_len[grep("_|M", names(chr_len), invert = TRUE)]
  my_ideo <- GenomicRanges::GRanges(
    seqnames = names(chr_len),
    ranges = IRanges::IRanges(start = 1, end = chr_len)
  )

  Seqinfo::seqlengths(my_ideo) <- width(my_ideo)
  Seqinfo::seqlevels(ranges, pruning.mode = "coarse") <-
    Seqinfo::seqlevels(my_ideo)
  Seqinfo::seqinfo(ranges) <- Seqinfo::seqinfo(my_ideo)

  g_hypo <- ranges[ranges$deltabetas < 0]
  g_hyper <- ranges[ranges$deltabetas > 0]

  if (length(g_hypo) > 0) values(g_hypo)$id <- "Down-regulated"
  if (length(g_hyper) > 0) values(g_hyper)$id <- "Up-regulated"

  p <- ggbio::ggbio() +
    ggbio::circle(
      my_ideo,
      geom = "ideo",
      fill = "gray92",
      color = "gray60",
      radius = 40,
      trackWidth = 4
    ) +
    ggbio::circle(
      my_ideo,
      geom = "text",
      aes(label = seqnames),
      vjust = 0.5,
      hjust = 0.5,
      radius = 40,
      trackWidth = 4,
      size = 2.2,
      color = "black",
      fontface = "bold"
    )

  if (length(g_hypo) > 0) {
    p <- p +
      ggbio::circle(
        g_hypo,
        geom = "point",
        size = point_size,
        aes(x = .data$midpoint, y = .data$deltabetas, color = .data$id),
        radius = 21,
        trackWidth = 18
      )
  }

  if (length(g_hyper) > 0) {
    p <- p +
      ggbio::circle(
        g_hyper,
        geom = "point",
        size = point_size,
        aes(x = .data$midpoint, y = .data$deltabetas, color = .data$id),
        radius = 45, trackWidth = 18
      )
  }

  p +
    ggplot2::scale_colour_manual(
      values = c("Up-regulated" = "#2E86AB", "Down-regulated" = "#A23B72"),
      name = "Direction"
    )
}
