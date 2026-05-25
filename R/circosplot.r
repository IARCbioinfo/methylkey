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
#' @export
#'
circosplot <- function(ranges, genome, label_probes = NULL,
                       label_size = 2.5, point_size = 1.2) {

  suppressPackageStartupMessages(require(ggbio))

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
      ranges$Probe_ID <- paste0("probe_", seq_len(length(ranges)))
    }
  }

  chr_len <- seqlengths(species)
  chr_len <- chr_len[grep("_|M", names(chr_len), invert = TRUE)]
  myIdeo <- GRanges(
    seqnames = names(chr_len),
    ranges = IRanges(start = 1, end = chr_len)
  )

  seqlengths(myIdeo) <- width(myIdeo)
  seqlevels(ranges, pruning.mode = "coarse") <- seqlevels(myIdeo)
  seqinfo(ranges) <- seqinfo(myIdeo)

  g.hypo <- ranges[ranges$deltabetas < 0]
  g.hyper <- ranges[ranges$deltabetas > 0]

  if (length(g.hypo) > 0) values(g.hypo)$id <- "Down-regulated"
  if (length(g.hyper) > 0) values(g.hyper)$id <- "Up-regulated"

  p <- ggbio() +
    circle(
      myIdeo,
      geom = "ideo",
      fill = "gray92",
      color = "gray60",
      radius = 40,
      trackWidth = 4
    ) +
    circle(
      myIdeo,
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

  if (length(g.hypo) > 0) {
    p <- p +
      circle(
        g.hypo,
        geom = "point",
        size = point_size,
        aes(x = midpoint, y = "deltabetas", color = id),
        radius = 21,
        trackWidth = 18
      )
  }

  if (length(g.hyper) > 0) {
    p <- p +
      circle(
        g.hyper,
        geom = "point",
        size = point_size,
        aes(x = midpoint, y = "deltabetas", color = id),
        radius = 45, trackWidth = 18
      )
  }

  p <- p +
    scale_colour_manual(
      values = c("Up-regulated" = "#2E86AB", "Down-regulated" = "#A23B72"),
      name = "Direction"
    )

  return(p)
}
