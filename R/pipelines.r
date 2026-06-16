#' Convert Illumina IDAT Files to Beta Values Using sesame
#'
#' This function takes Illumina IDAT files,
#'   processes them using the sesame package,
#'   and returns beta values (DNA methylation values)
#'   along with quality control statistics and plots.
#'
#' @param idat The path to the directory containing Illumina IDAT files.
#' @param prep The sesame preprocessing method to use. Default is "QCDPB".
#' @param sampleSheet A data frame containing sample information,
#'   including barcode, sex, and age.
#' @param na The maximum proportion of missing values allowed in beta values.
#'   Must be between 0 and 1. Default is 0.2.
#' @param ncore The number of CPU cores to use for parallel processing.
#'   Default is 4.
#'
#' @return An object of class "MethyLumiMethyLumiSet" containing beta values
#'   and associated metadata, including quality control statistics and plots.
#'
#' @importFrom assertthat assert_that
#' @importFrom sesame openSesame inferSex getBetas
#' @importFrom sesame sesameQC_calcStats sesameQC_getStats
#' @importFrom BiocParallel MulticoreParam
#' @importFrom dplyr left_join
#' @importFrom tibble tibble
#'
#' @export
sesame2betas <- function(
  idat = NULL,
  prep = "QCDPB",
  sample_sheet = NULL,
  na = 0.2,
  ncore = 4,
  clock_models = NULL
) {

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' is required for this function.",
      "Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("sesameData", quietly = TRUE)) {
    stop("Package 'sesameData' is required for this function.",
      "Please install it.",
      call. = FALSE
    )
  }

  assertthat::assert_that(
    dir.exists(idat),
    msg = paste0("Directory not found: ", idat)
  )
  assertthat::assert_that(
    is.data.frame(sample_sheet),
    msg = "sample_sheet must be a data frame"
  )
  assertthat::assert_that(
    na >= 0 && na <= 1,
    msg = "na must be between 0 and 1"
  )

  # Format sampleSheet
  sample_sheet <- format_sample_sheet(sample_sheet)

  # Load IDAT into sdfs list
  sdfs <- sesame::openSesame(
    idat,
    prep = prep,
    func = NULL,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )

  # Save sdfs object to avoid recalculating if needed later
  saveRDS(sdfs, file = file.path("sdfs.rds"))

  # Extract betas matrix
  betas <- sesame::openSesame(
    sdfs,
    prep = "",
    func = sesame::getBetas,
    mask = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )

  # Infer sex if betas are available
  if (ncol(betas) > 0) {
    infered_sex <- sapply(colnames(betas), function(barcode) {
      sesame::inferSex(betas[, barcode])
    })
    infered_sex <- tibble::tibble(
      barcode = names(infered_sex),
      infered_sex = infered_sex
    )
    sample_sheet <- dplyr::left_join(
      sample_sheet,
      infered_sex,
      by = "barcode"
    )
  }

  # Infer age using clock models if provided
  if (!is.null(clock_models) && length(clock_models) > 0) {
    for (model_name in names(clock_models)) {
      model_file <- clock_models[[model_name]]
      if (file.exists(model_file)) {
        model <- readr::read_rds(model_file)
        na_fallback <- "na_fallback" %in% names(model$param)
        sample_sheet[[model_name]] <- sapply(
          sample_sheet$barcode, function(barcode) {
            sesame::predictAge(betas[, barcode], model, na_fallback)
          }
        )
      }
    }
  }

  # Create methylkey Beta object
  meth <- new_betas(betas, sample_sheet, na)
  metadata(meth)$betas.pipeline.name <- "sesame"
  metadata(meth)$preprocessing.method <- prep

  # Add QC statistics
  qcs <- sesame::openSesame(
    sdfs,
    prep = "",
    func = sesame::sesameQC_calcStats,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )
  metadata(meth)$qcs <- purrr::map(qcs, ~sesame::sesameQC_getStats(.)) |>
    dplyr::bind_rows(.id = "barcode")

  meth
}

#' Convert Illumina IDAT Files to Beta Values Using minfi
#'
#' This function takes Illumina IDAT files,
#'   processes them using the minfi package,
#'   and returns beta values (DNA methylation values)
#'   along with quality control statistics and plots.
#'
#' @param idat The path to the directory containing Illumina IDAT files.
#' @param sample_sheet A data frame containing sample information,
#'   including Basename, sex, and age.
#' @param na The maximum proportion of missing values allowed in beta values.
#'   Must be between 0 and 1. Default is 0.2.
#' @param compositeCellType A composite cell type to estimate cell counts,
#'   e.g., "Blood," "CordBloodCombined," etc.
#'
#' @return An object of class "MethyLumiMethyLumiSet" containing beta values
#'   and associated metadata, including quality control statistics and plots.
#'
#' @importFrom assertthat assert_that
#'
#' @export
minfi2betas <- function(
  idat = NULL,
  sample_sheet = NULL,
  na = 0.2,
  composite_cell_type = "",
  pval = 0.2
) {

  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package 'minfi' is required for this function.",
      "Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("wateRmelon", quietly = TRUE)) {
    stop("Package 'wateRmelon' is required for this function.",
      "Please install it.",
      call. = FALSE
    )
  }

  composite_cell_type_valid <- c(
    "Blood", "CordBloodCombined", "BloodExtended", "CordBlood",
    "CordBloodNorway", "CordTissueAndBlood", "DLPFC"
  )
  cell_types_standard <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
  #cell_types_extended <- c(
  #  "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv",
  #  "Eos", "Mono", "Neu", "NK", "Treg"
  #)

  assertthat::assert_that(
    dir.exists(idat),
    msg = paste0("Directory not found: ", idat)
  )
  assertthat::assert_that(
    "Basename" %in% colnames(sample_sheet),
    msg = "Basename column required in sample_sheet"
  )
  assertthat::assert_that(
    is.data.frame(sample_sheet),
    msg = "sample_sheet must be a data frame"
  )
  assertthat::assert_that(
    na >= 0 && na <= 1,
    msg = "na must be between 0 and 1"
  )
  assertthat::assert_that(
    pval >= 0 && pval <= 1,
    msg = "pval must be between 0 and 1"
  )

  rgset <- minfi::read.metharray.exp(
    base = idat,
    targets = sample_sheet,
    force = TRUE
  )

  gmsetex <- minfi::mapToGenome(rgset)
  estsex <- minfi::getSex(gmsetex)
  sample_sheet$infered_sex <- estsex$predictedSex

  mset <- minfi::preprocessRaw(rgset)
  mset <- minfi::fixMethOutliers(mset)
  qcs <- minfi::getQC(mset)

  # Normalization
  betas <- minfi::getBeta(rgset)
  isna <- is.na(betas)
  mset <- minfi::preprocessFunnorm(rgset, sex = estsex$predictedSex)
  betas <- minfi::getBeta(mset)

  # Restore NA status after normalization
  isna <- isna[match(rownames(betas), rownames(isna)), ]
  betas[which(isna)] <- NA

  pvalues <- minfi::detectionP(rgset)
  betas[pvalues[rownames(betas), ] > pval] <- NA

  # Estimate cell counts if requested
  if (composite_cell_type != "" &&
      composite_cell_type %in% composite_cell_type_valid
  ) {

    if (get_plateform(betas) == "IlluminaHumanMethylationEPIC") {
      if (!requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE)) {
        warning("FlowSorted.Blood.EPIC package",
                "required for cell count estimation. Skipping.")
      } else {
        cc <- FlowSorted.Blood.EPIC::estimateCellCounts2(
          rgset,
          composite_cell_type,
          cell_types = cell_types_standard,
          referencePlatform = get_plateform(betas)
        )
        sample_sheet <- cbind(sample_sheet, cc$prop)
        sample_sheet$nlr <- sample_sheet$neu / sample_sheet$nk
      }
    } else if (get_plateform(betas) %in% c(
      "IlluminaHumanMethylation450k",
      "IlluminaHumanMethylation27k"
    )) {
      cc <- minfi::estimateCellCounts(
        rgset,
        composite_cell_type,
        cell_types = cell_types_standard,
        referencePlatform = get_plateform(betas)
      )
      sample_sheet <- cbind(sample_sheet, cc)
      sample_sheet$nlr <- sample_sheet$Neu / sample_sheet$NK
    } else {
      warning(sprintf(
        "Cell count estimation not available for platform: %s",
        get_plateform(betas)
      ))
    }

    # Define immune score
    if ("nlr" %in% colnames(sample_sheet)) {
      breaks <- c(-Inf, 0.7, 1, 2, 3, Inf)
      labels <- c("bad", "average", "good", "average", "bad")
      sample_sheet$epimmune <- cut(sample_sheet$nlr,
        breaks = breaks,
        labels = labels
      )
    }
  }

  sample_sheet <- cbind(
    sample_sheet,
    wateRmelon::agep(betas, method = "all")
  )

  meth <- new_betas(betas, sample_sheet, na)
  metadata(meth)$betas.pipeline.name <- "minfi"
  metadata(meth)$qcs <- qcs
  metadata(meth)$platform <- get_plateform(betas)
  metadata(meth)$celltype_estimation <- composite_cell_type

  meth
}