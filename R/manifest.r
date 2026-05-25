
getZhouManifest <- function(plateform) {

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
}

getZhouAnnotation <- function(plateform) {

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

  annotation <- annotation |> dplyr::rename(Probe_ID = probeID)
}

getIlluminaAnnotation <- function(plateform) {

  switch(
    plateform,

    "IlluminaHumanMethylation450k" = {
      cbind(
        IlluminaHumanMethylation450kanno.ilmn12.hg19::Islands.UCSC,
        IlluminaHumanMethylation450kanno.ilmn12.hg19::Other
      ) |> as.data.frame() |>
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
      ) |> as.data.frame() |>
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
      ) |> as.data.frame() |>
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
      ) |> as.data.frame() |>
        dplyr::select("GeneName_UCSC",
                      "Transcript_UCSC",
                      "Feature_UCSC",
                      "Islands_Name",
                      "Relation_to_Island")
    }

  ) |>
    tibble::rownames_to_column("Probe_ID")

}


#' getAnnotedManifest
#'
#' retrieve infinium annotations
#'
#' @param plateform plateform
#'
#' @return annotated manifest
#'
#' @export
#'
getManifest <- function(plateform) {

  manifest <- getZhouManifest(plateform) |>
    dplyr::rename(chr = CpG_chrm, pos =  CpG_beg, strand = mapYD_A) |>
    dplyr::select(Probe_ID, chr, pos, strand) |>
    mutate(pos = pos + 1) |>
    mutate(strand = ifelse(is.na(strand), NA, ifelse(strand == "f", "+", "-")))

  return(manifest)
}


get_annotated_manifest <- function(se, annot = NULL) {

  pl <- methM@metadata$plateform

  if (is.null(annot)) {
    annot1 = getZhouAnnotation(pl)
    annot2 = getIlluminaAnnotation(pl)
    annot <- dplyr::full_join(annot1, annot2, by = "Probe_ID")
  }

  annotated_manifest <-
    annot[match(rowData(se)$Probe_ID, annot$Probe_ID), ] |>
    dplyr::filter(!is.na(Probe_ID)) |>
    S4Vectors::DataFrame()

}