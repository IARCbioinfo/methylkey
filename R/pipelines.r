#' Convert Illumina IDAT Files to Beta Values Using sesame
#'
#' This function takes Illumina IDAT files, processes them using the sesame package, and returns beta values (DNA methylation values) along with quality control statistics and plots.
#'
#' @param idat The path to the directory containing Illumina IDAT files.
#' @param prep The sesame preprocessing method to use. Default is "QCDPB".
#' @param sampleSheet A data frame containing sample information, including barcode, sex, and age.
#' @param na The maximum proportion of missing values allowed in beta values. Must be between 0 and 1. Default is 0.2.
#' @param ncore The number of CPU cores to use for parallel processing. Default is 4.
#'
#' @return An object of class "MethyLumiMethyLumiSet" containing beta values and associated metadata, including quality control statistics and plots.
#'
#' @export
sesame2Betas <- function(
  idat = NULL,
  prep = "QCDPB",
  sampleSheet = NULL,
  na = 0.2,
  ncore = 4,
  Clock_models = NULL
) {

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' is required for this function. Please install it.")
  }
  if (!requireNamespace("sesameData", quietly = TRUE)) {
    stop("Package 'sesameData' is required for this function. Please install it.")
  }

  assertthat::assert_that(
    dir.exists(idat),
    msg = paste0("Directory not found: ", idat)
  )
  assertthat::assert_that(
    is.data.frame(sampleSheet),
    msg = "sampleSheet must be a data frame"
  )
  assertthat::assert_that(na >= 0 && na <= 1, msg = "na must be between 0 and 1")

  # Format sampleSheet
  sampleSheet <- formatSampleSheet(sampleSheet)

  # Load IDAT into sdfs list
  message("Loading IDAT files from sesame...")
  sdfs <- sesame::openSesame(
    idat,
    prep = prep,
    func = NULL,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )

  # Save sdfs object to avoid recalculating if needed later
  saveRDS(sdfs, file = file.path("sdfs.rds"))

  # Extract betas matrix
  message("Extracting beta values...")
  betas <- sesame::openSesame(
    sdfs,
    prep = "",
    func = sesame::getBetas,
    mask = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )

  # Infer sex if betas are available
  if (ncol(betas) > 0) {
    message("Inferring sample sex...")
    inferedSex <- sapply(colnames(betas), function(barcode) {
      sesame::inferSex(betas[, barcode])
    })
    inferedSex <- tibble::tibble(
      barcode = names(inferedSex),
      inferedSex = inferedSex
    )
    sampleSheet <- dplyr::left_join(sampleSheet, inferedSex, by = "barcode")
  }

  # Infer age using clock models if provided
  if (!is.null(Clock_models) && length(Clock_models) > 0) {
    message("Computing epigenetic age...")
    for (model_name in names(Clock_models)) {
      model_file <- Clock_models[[model_name]]
      if (file.exists(model_file)) {
        model <- readr::read_rds(model_file)
        na_fallback <- "na_fallback" %in% names(model$param)
        sampleSheet[[model_name]] <- sapply(sampleSheet$barcode, function(barcode) {
          sesame::predictAge(betas[, barcode], model, na_fallback)
        })
      }
    }
  }

  # Create methylkey Beta object
  message("Creating methylkey object...")
  meth <- newBetas(betas, sampleSheet, na)
  metadata(meth)$betas.pipeline.name <- "sesame"
  metadata(meth)$preprocessing.method <- prep

  # Add QC statistics
  message("Computing QC statistics...")
  qcs <- sesame::openSesame(
    sdfs,
    prep = "",
    func = sesame::sesameQC_calcStats,
    BPPARAM = BiocParallel::MulticoreParam(ncore)
  )
  metadata(meth)$qcs <- purrr::map(qcs, ~sesame::sesameQC_getStats(.)) |>
    dplyr::bind_rows(.id = "barcode")

  message("Sesame processing complete.")
  return(meth)
}

#' Convert Illumina IDAT Files to Beta Values Using minfi
#'
#' This function takes Illumina IDAT files, processes them using the minfi package, and returns beta values (DNA methylation values) along with quality control statistics and plots.
#'
#' @param idat The path to the directory containing Illumina IDAT files.
#' @param sampleSheet A data frame containing sample information, including Basename, sex, and age.
#' @param na The maximum proportion of missing values allowed in beta values. Must be between 0 and 1. Default is 0.2.
#' @param compositeCellType A composite cell type to estimate cell counts, e.g., "Blood," "CordBloodCombined," etc.
#'
#' @return An object of class "MethyLumiMethyLumiSet" containing beta values and associated metadata, including quality control statistics and plots.
#'
#' @export
minfi2Betas <- function(
  idat = NULL,
  sampleSheet = NULL,
  na = 0.2,
  compositeCellType = "",
  pval = 0.2
) {

  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package 'minfi' is required for this function. Please install it.")
  }
  if (!requireNamespace("wateRmelon", quietly = TRUE)) {
    stop("Package 'wateRmelon' is required for this function. Please install it.")
  }

  compositeCellType_valid <- c(
    "Blood", "CordBloodCombined", "BloodExtended", "CordBlood",
    "CordBloodNorway", "CordTissueAndBlood", "DLPFC"
  )
  cellTypes_standard <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
  cellTypes_extended <- c(
    "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv",
    "Eos", "Mono", "Neu", "NK", "Treg"
  )

  assertthat::assert_that(
    dir.exists(idat),
    msg = paste0("Directory not found: ", idat)
  )
  assertthat::assert_that(
    "Basename" %in% colnames(sampleSheet),
    msg = "Basename column required in sampleSheet"
  )
  assertthat::assert_that(
    is.data.frame(sampleSheet),
    msg = "sampleSheet must be a data frame"
  )
  assertthat::assert_that(na >= 0 && na <= 1, msg = "na must be between 0 and 1")
  assertthat::assert_that(pval >= 0 && pval <= 1, msg = "pval must be between 0 and 1")

  message("Loading IDAT files using minfi...")
  RGset <- minfi::read.metharray.exp(
    base = idat,
    targets = sampleSheet,
    force = TRUE
  )

  message("Mapping to genome...")
  GMsetEx <- minfi::mapToGenome(RGset)
  estSex <- minfi::getSex(GMsetEx)
  sampleSheet$inferedsex <- estSex$predictedSex

  message("Preprocessing...")
  MSet <- minfi::preprocessRaw(RGset)
  MSet <- minfi::fixMethOutliers(MSet)
  qcs <- minfi::getQC(MSet)

  # Normalization
  message("Normalizing...")
  betas <- minfi::getBeta(RGset)
  isna <- is.na(betas)
  MSet <- minfi::preprocessFunnorm(RGset, sex = estSex$predictedSex)
  betas <- minfi::getBeta(MSet)

  # Restore NA status after normalization
  isna <- isna[match(rownames(betas), rownames(isna)), ]
  betas[which(isna)] <- NA

  message("Applying detection p-value filter...")
  pvalues <- minfi::detectionP(RGset)
  n_filtered <- sum(pvalues > pval)
  betas[pvalues[rownames(betas), ] > pval] <- NA
  message(
    sprintf(
      "Marked %d low quality probes as NA (p-value > %.3f)",
      n_filtered,
      pval
    )
  )

  # Estimate cell counts if requested
  if (compositeCellType != "" && compositeCellType %in% compositeCellType_valid) {
    message(sprintf("Estimating cell counts for: %s", compositeCellType))

    if (getPlateform(betas) == "IlluminaHumanMethylationEPIC") {
      if (!requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE)) {
        warning("FlowSorted.Blood.EPIC package required for cell count estimation. Skipping.")
      } else {
        cc <- minfi::estimateCellCounts2(RGset, compositeCellType,
          cellTypes = cellTypes_standard,
          referencePlatform = getPlateform(betas)
        )
        sampleSheet <- cbind(sampleSheet, cc$prop)
        sampleSheet$nlr <- sampleSheet$neu / sampleSheet$nk
      }
    } else if (getPlateform(betas) %in% c(
      "IlluminaHumanMethylation450k",
      "IlluminaHumanMethylation27k"
    )) {
      cc <- minfi::estimateCellCounts(RGset, compositeCellType,
        cellTypes = cellTypes_standard,
        referencePlatform = getPlateform(betas)
      )
      sampleSheet <- cbind(sampleSheet, cc)
      sampleSheet$nlr <- sampleSheet$Neu / sampleSheet$NK
    } else {
      warning(sprintf(
        "Cell count estimation not available for platform: %s",
        getPlateform(betas)
      ))
    }

    # Define immune score
    if ("nlr" %in% colnames(sampleSheet)) {
      breaks <- c(-Inf, 0.7, 1, 2, 3, Inf)
      labels <- c("bad", "average", "good", "average", "bad")
      sampleSheet$epimmune <- cut(sampleSheet$nlr,
        breaks = breaks,
        labels = labels
      )
    }
  }

  message("Computing epigenetic age predictors...")
  sampleSheet <- cbind(sampleSheet, wateRmelon::agep(betas, method = "all"))

  message("Creating methylkey object...")
  meth <- newBetas(betas, sampleSheet, na)
  metadata(meth)$betas.pipeline.name <- "minfi"
  metadata(meth)$qcs <- qcs
  metadata(meth)$platform <- getPlateform(betas)
  metadata(meth)$celltype_estimation <- compositeCellType

  message("Minfi processing complete.")
  return(meth)
}




