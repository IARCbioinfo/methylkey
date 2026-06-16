
#' getZhouManifest
#'
#' @param plateform Type of plateform
#' @return A data.frame manifest
#' @importFrom readr read_tsv
get_zhou_manifest <- function(plateform) {

  manifest <- NULL

  zhou_lab_url <-
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/"

  switch(plateform,
    "IlluminaHumanMethylation450k" = {
      manifest <- readr::read_tsv(
        paste0(zhou_lab_url, "HM450/HM450.hg38.manifest.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPIC" = {
      manifest <- readr::read_tsv(
        paste0(zhou_lab_url, "EPIC/EPIC.hg38.manifest.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPIC+" = {
      manifest <- readr::read_tsv(
        paste0(zhou_lab_url, "EPIC+/EPIC+.hg38.manifest.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPICv2" = {
      manifest <- readr::read_tsv(
        paste0(zhou_lab_url, "EPICv2/EPICv2.hg38.manifest.tsv.gz")
      )
    },
    "IlluminaMouseMethylation285k" = {
      manifest <- readr::read_tsv(
        paste0(zhou_lab_url, "MM285/MM285.mm10.manifest.tsv.gz")
      )
    }
  )
  manifest
}


#' get_zhou_annotation
#'
#' @param plateform Type of plateform
#' @return A data.frame annotation
#' @importFrom readr read_tsv
#' @importFrom dplyr rename
get_zhou_annotation <- function(plateform) {

  annotation <- NULL

  zhou_lab_url <-
    "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/"

  switch(plateform,
    "IlluminaHumanMethylation450k" = {
      annotation <- readr::read_tsv(
        paste0(zhou_lab_url, "HM450/HM450.hg38.manifest.gencode.v37.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPIC" = {
      annotation <- readr::read_tsv(
        paste0(zhou_lab_url, "EPIC/EPIC.hg38.manifest.gencode.v37.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPIC+" = {
      annotation <- readr::read_tsv(
        paste0(zhou_lab_url, "EPIC%2B/EPIC%2B.hg38.manifest.gencode.v41.tsv.gz")
      )
    },
    "IlluminaHumanMethylationEPICv2" = {
      annotation <- readr::read_tsv(
        paste0(zhou_lab_url, "EPICv2/EPICv2.hg38.manifest.gencode.v41.tsv.gz")
      )
    },
    "IlluminaMouseMethylation285k" = {
      annotation <- readr::read_tsv(
        paste0(zhou_lab_url, "MM285/MM285.mm10.manifest.gencode.vM25.tsv.gz")
      )
    }
  )

  annotation |> dplyr::rename(Probe_ID = "probeID")

}


#' get_illumina_annotation
#'
#' @param plateform Type of plateform
#' @return A data.frame annotation
#' @importFrom dplyr select mutate
#' @importFrom tibble rownames_to_column
get_illumina_annotation <- function(plateform) {

  switch(
    plateform,

    "IlluminaHumanMethylation450k" = {
      cbind(
        IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC,
        IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
      ) |>
        as.data.frame() |>
        dplyr::select("UCSC_RefGene_Name",
                      "UCSC_RefGene_Accession",
                      "UCSC_RefGene_Group",
                      "Islands_Name",
                      "Relation_to_Island")
    },

    "IlluminaHumanMethylationEPIC" = {
      cbind(
        IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Islands.UCSC,
        IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other
      ) |>
        as.data.frame() |>
        dplyr::select("UCSC_RefGene_Name",
                      "UCSC_RefGene_Accession",
                      "UCSC_RefGene_Group",
                      "Islands_Name",
                      "Relation_to_Island")
    },

    "IlluminaHumanMethylationEPICv2" = {
      cbind(
        IlluminaHumanMethylationEPICv2anno.20a1.hg38::Islands.UCSC,
        IlluminaHumanMethylationEPICv2anno.20a1.hg38::Other
      ) |>
        as.data.frame() |>
        dplyr::select("UCSC_RefGene_Name",
                      "UCSC_RefGene_Accession",
                      "UCSC_RefGene_Group",
                      "Islands_Name",
                      "Relation_to_Island")
    },

    "IlluminaMouseMethylation285k" = {
      cbind(
        IlluminaMouseMethylationanno.12.v1.mm10::Islands.UCSC,
        IlluminaMouseMethylationanno.12.v1.mm10::GenesUCSC
      ) |>
        as.data.frame() |>
        dplyr::select("GeneName_UCSC",
                      "Transcript_UCSC",
                      "Feature_UCSC",
                      "Islands_Name",
                      "Relation_to_Island")
    }

  ) |>
    tibble::rownames_to_column("Probe_ID")

}


#' getManifest
#'
#' retrieve infinium annotations
#'
#' @param plateform plateform
#' @return annotated manifest
#' @importFrom dplyr rename select mutate
#' @export
get_manifest <- function(plateform) {

  get_zhou_manifest(plateform) |>
    dplyr::rename(chr = "CpG_chrm", pos = "CpG_beg", strand = "mapYD_A") |>
    dplyr::select(.data$Probe_ID, .data$chr, .data$pos, .data$strand) |>
    dplyr::mutate(pos = .data$pos + 1) |>
    dplyr::mutate(
      strand = ifelse(
        is.na(.data$strand),
        NA,
        ifelse(.data$strand == "f", "+", "-")
      )
    )
}


#' get_annotated_manifest
#'
#' @param se SummarizedExperiment
#' @param annot annotation optionnelle
#' @return DataFrame annoté
#' @importFrom dplyr full_join filter
#' @importFrom S4Vectors DataFrame
get_annotated_manifest <- function(se, annot = NULL) {

  pl <- se@metadata$plateform

  if (is.null(annot)) {
    annot1 = get_zhou_annotation(pl)
    annot2 = get_illumina_annotation(pl)
    annot <- dplyr::full_join(annot1, annot2, by = "Probe_ID")
  }

  annot[match(SummarizedExperiment::rowData(se)$Probe_ID, annot$Probe_ID), ] |>
    dplyr::filter(!is.na(.data$Probe_ID)) |>
    S4Vectors::DataFrame()

}