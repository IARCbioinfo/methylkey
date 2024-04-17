#' Class Betas
#'
#' This class extend SummarizedExperiment, with new methods, it made it specific for betas values.
#'
#' @seealso
#' \code{\link{SummarizedExperiment-class}} : The parent class
#'
#' @family data container
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
setClass("Betas",
         slots=list(),
         contains="SummarizedExperiment",
)

#' Class Betas
#'
#' This class extend SummarizedExperiment, with new methods, it made it specific for mvalues 
#'
#' @seealso
#' \code{\link{SummarizedExperiment-class}} : The parent class
#'
#' @family data container
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' 
#' @export
setClass("Mvals",
         slots=list(),
         contains="SummarizedExperiment",
)

#' Generic function to extract column data from an object
#'
#' @param x An object from which to extract column data.
#'
#' @return The extracted column data.
#'
#' @export
setGeneric("colData", function(x) standardGeneric("colData") )

#' Method for extracting column data from a "Betas" object
#'
#' This method extracts column data from a "Betas" object and returns it as a data frame.
#'
#' @param x A "Betas" object.
#'
#' @return A data frame containing the extracted column data.
#'
#' @export
setMethod("colData", signature(x="Betas"),definition=function(x) data.frame(x@colData))

#' Method for extracting column data from a "Mvals" object
#'
#' This method extracts column data from a "Mvals" object and returns it as a data frame.
#'
#' @param x A "Mvals" object.
#'
#' @return A data frame containing the extracted column data.
#'
#' @export
setMethod("colData", signature(x="Mvals"),definition=function(x) data.frame(x@colData))

#' Method for extracting column data from a "NULL" object
#'
#' This method is used to extract column data from a "NULL" object. In this case, it always returns NULL, indicating that there is no column data to extract from a NULL object.
#'
#' @param x A "NULL" object.
#'
#' @return NULL, as there is no column data to extract from a "NULL" object.
#'
#' @export
setMethod("colData", signature(x="NULL"), definition = function(x) NULL)

#' Generic function to set column data for an object
#'
#' This is a generic function to set column data for an object. It defines
#' the interface for setting column data and can be used with different classes.
#'
#' @param x An object to which to set the column data.
#' @param pdata A data frame containing the new column data.
#'
#' @return The modified object with the updated column data.
#'
#' @export
setGeneric("setColData", function(x, pdata) standardGeneric("setColData") )

#' Method for setting column data for a "Betas" object
#'
#' This method sets the column data for a "Betas" object using the provided data frame.
# It modifies the object in place and returns the modified object.
#'
#' @param x A "Betas" object to which to set the column data.
#' @param pdata A data frame containing the new column data.
#'
#' @return The modified "Betas" object with the updated column data.
#'
#' @export
setMethod("setColData", signature("Betas"), function(x, pdata){
  x@colData <- S4Vectors::DataFrame(pdata)
  x
})

#' Generic function to retrieve beta values from a "Betas" object
#'
#' This is a generic function to retrieve beta values from a "Betas" object.
# It provides options to filter the data based on certain conditions.
#'
#' @param x A "Betas" object from which to retrieve beta values.
#' @param masked Logical. If FALSE (default), filter out masked values.
#' @param na Logical. If FALSE (default), filter out NA values.
#' @param sex Logical. If TRUE (default), retain sex chromosomes data.
#'
#' @return A matrix of beta values, optionally filtered based on the specified conditions.
#'
#' @export
setGeneric("getBetas", function(x,masked=FALSE, na=FALSE, sex=TRUE) standardGeneric("getBetas"))

#' Method for retrieving beta values from a "Betas" object
#'
#' This method retrieves beta values from a "Betas" object, optionally filtering the data based on the specified conditions.
#'
#' @param x A "Betas" object from which to retrieve beta values.
#' @param masked Logical. If FALSE (default), filter out masked values.
#' @param na Logical. If FALSE (default), filter out NA values.
#' @param sex Logical. If TRUE (default), retain sex chromosomes data.
#'
#' @return A matrix of beta values, optionally filtered based on the specified conditions.
#'
#' @export
setMethod("getBetas", signature("Betas"),
  definition = function(x, masked = FALSE, na=FALSE, sex = TRUE) {
    if(!masked & !na & sex) return(assays(x)$betas)
    mask_general_condition <- if (!masked) { rowData(x)$M_general == FALSE } else { TRUE }
    na_condition <- if (!na) { rowData(x)$Na == FALSE } else { TRUE }
    mask_chrm_condition <- if (!sex) { rowData(x)$MASK_chrm == FALSE } else { TRUE }
    filtered_data <- assays(x)$betas[which(mask_general_condition & na_condition & mask_chrm_condition), ]
    return(filtered_data)
  }
)

#' Create a generic function for getting beta value ranges
#'
#' This function sets up a generic function, `getBetasRanges`, for getting beta value ranges from a Betas object.
#'
#' @title Create generic function for getting beta value ranges
#' @param x A Betas object.
#' @param masked Logical indicating whether to include masked CpGs (default is FALSE).
#' @param na Logical indicating whether to include CpGs marked as NA (default is FALSE).
#' @param sex Logical indicating whether to include sex chromosomes (default is TRUE).
#' 
#' @export
#' 
setGeneric("getBetasRanges", function(x,masked=FALSE, na=FALSE, sex=TRUE) standardGeneric("getBetasRanges"))

setMethod("getBetasRanges", signature("Betas"),
          definition = function(x, masked = FALSE, na=FALSE, sex = TRUE) {
            if(!masked & !na & sex) return(assays(x)$betas)
            mask_general_condition <- if (!masked) { rowData(x)$M_general == FALSE } else { TRUE }
            na_condition <- if (!na) { rowData(x)$Na == FALSE } else { TRUE }
            mask_chrm_condition <- if (!sex) { rowData(x)$MASK_chrm == FALSE } else { TRUE }
            df = bind_cols(as.data.frame(rowData(x)), assays(x)$betas) %>% mutate()
            filtered_data <- df[which(mask_general_condition & na_condition & mask_chrm_condition), ]
            filtered_data <- GenomicRanges::makeGRangesFromDataFrame(filtered_data, seqnames.field="CpG_chrm", start.field="CpG_beg", end.field="CpG_end", keep.extra.columns = TRUE )
            return(filtered_data)
          }
)


#' Generic function to retrieve M-values from an "Mvals" object
#'
#' This is a generic function to retrieve M-values from an "Mvals" object.
# It provides a consistent interface for accessing M-values across different classes.
#'
#' @param x An "Mvals" object from which to retrieve M-values.
#'
#' @return A matrix of M-values.
#'
#' @export
setGeneric("getMvals", function(x, ...) standardGeneric("getMvals"))

#' Method for retrieving M-values from an "Mvals" object
#'
#' This method retrieves M-values from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve M-values.
#'
#' @return A matrix of M-values.
#'
#' @export
setMethod("getMvals", signature("Mvals"), definition = function(x) { return(assays(x)$mvals) })

#' Generic function to retrieve sample information from an object
#'
#' This is a generic function to retrieve sample information from an object.
# It provides a consistent interface for accessing sample information across different classes.
#'
#' @param x An object from which to retrieve sample information.
#'
#' @return Sample information, typically as a data frame or other suitable format.
#'
#' @export
setGeneric("samples", function(x) standardGeneric("samples"))

#' Method for retrieving sample information from a "Betas" object
#'
#' This method retrieves sample information from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve sample information.
#'
#' @return Sample information as a data frame.
#'
#' @export
setMethod("samples",signature("Betas"),definition=function(x) x@colData$samples )

#' Method for retrieving sample information from an "Mvals" object
#'
#' This method retrieves sample information from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve sample information.
#'
#' @return Sample information as a data frame.
#'
#' @export
setMethod("samples",signature("Mvals"),definition=function(x) x@colData$samples )


#' Generic function to retrieve the number of samples from an object
#'
#' This is a generic function to retrieve the number of samples from an object.
# It provides a consistent interface for accessing the number of samples across different classes.
#'
#' @param x An object from which to retrieve the number of samples.
#'
#' @return The number of samples.
#'
#' @export
setGeneric("nbsamples", function(x) standardGeneric("nbsamples"))

#' Method for retrieving the number of samples from a "Betas" object
#'
#' This method retrieves the number of samples from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve the number of samples.
#'
#' @return The number of samples in the "Betas" object.
#'
#' @export
setMethod("nbsamples",signature("Betas"),definition=function(x) nrow(x@colData) )

#' Method for retrieving the number of samples from an "Mvals" object
#'
#' This method retrieves the number of samples from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve the number of samples.
#'
#' @return The number of samples in the "Mvals" object.
#'
#' @export
setMethod("nbsamples",signature("Mvals"),definition=function(x) nrow(x@colData) )

#' Generic function to retrieve variable names from an object
#'
#' This is a generic function to retrieve variable names (column names) from an object.
# It provides a consistent interface for accessing variable names across different classes.
#'
#' @param x An object from which to retrieve variable names.
#'
#' @return A character vector containing variable names.
#'
#' @export
setGeneric("variables", function(x) standardGeneric("variables"))

#' Method for retrieving variable names from a "Betas" object
#'
#' This method retrieves variable names (column names) from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve variable names.
#'
#' @return A character vector containing variable names.
#'
#' @export
setMethod("variables",signature("Betas"),definition=function(x) colnames(x@colData) )

#' Method for retrieving variable names from an "Mvals" object
#'
#' This method retrieves variable names (column names) from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve variable names.
#'
#' @return A character vector containing variable names.
#'
#' @export
setMethod("variables",signature("Mvals"),definition=function(x) colnames(x@colData) )

#' Generic function to retrieve Sentrix IDs from an object
#'
#' This is a generic function to retrieve Sentrix IDs from an object.
# It provides a consistent interface for accessing Sentrix IDs across different classes.
#'
#' @param x An object from which to retrieve Sentrix IDs.
#'
#' @return A character vector containing Sentrix IDs.
#'
#' @export
setGeneric("sentrix", function(x) standardGeneric("sentrix"))

#' Method for retrieving Sentrix IDs from a "Betas" object
#'
#' This method retrieves Sentrix IDs from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve Sentrix IDs.
#'
#' @return A character vector containing Sentrix IDs.
#'
#' @export
setMethod("sentrix",signature("Betas"),definition=function(x) colnames(assays(x)$betas) )

#' Method for retrieving Sentrix IDs from an "Mvals" object
#'
#' This method retrieves Sentrix IDs from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve Sentrix IDs.
#'
#' @return A character vector containing Sentrix IDs.
#'
#' @export
setMethod("sentrix",signature("Mvals"),definition=function(x) colnames(assays(x)$mvals) )

#' Generic function to retrieve the number of Sentrix IDs from an object
#'
#' This is a generic function to retrieve the number of Sentrix IDs from an object.
# It provides a consistent interface for accessing the number of Sentrix IDs across different classes.
#'
#' @param x An object from which to retrieve the number of Sentrix IDs.
#'
#' @return The number of Sentrix IDs as an integer.
#'
#' @export
setGeneric("nbsentrix", function(x) standardGeneric("nbsentrix"))

#' Method for retrieving the number of Sentrix IDs from a "Betas" object
#'
#' This method retrieves the number of Sentrix IDs from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve the number of Sentrix IDs.
#'
#' @return The number of Sentrix IDs as an integer.
#'
#' @export
setMethod("nbsentrix",signature("Betas"),definition=function(x) ncol(assays(x)$betas) )

#' Method for retrieving the number of Sentrix IDs from an "Mvals" object
#'
#' This method retrieves the number of Sentrix IDs from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve the number of Sentrix IDs.
#'
#' @return The number of Sentrix IDs as an integer.
#'
#' @export
setMethod("nbsentrix",signature("Mvals"),definition=function(x) ncol(assays(x)$mvals) )

#' Generic function to retrieve probe names from an object
#'
#' This is a generic function to retrieve probe names from an object.
# It provides a consistent interface for accessing probe names across different classes.
#'
#' @param x An object from which to retrieve probe names.
#'
#' @return A character vector containing probe names.
#'
#' @export
setGeneric("probes", function(x) standardGeneric("probes"))

#' Method for retrieving probe names from a "Betas" object
#'
#' This method retrieves probe names from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve probe names.
#'
#' @return A character vector containing probe names.
#'
#' @export
setMethod("probes",signature("Betas"),definition=function(x) rownames(assays(x)$betas) )

#' Method for retrieving probe names from an "Mvals" object
#'
#' This method retrieves probe names from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve probe names.
#'
#' @return A character vector containing probe names.
#'
#' @export
setMethod("probes",signature("Mvals"),definition=function(x) rownames(assays(x)$mvals) )

#' Generic function to retrieve the number of probes from an object
#'
#' This is a generic function to retrieve the number of probes from an object.
# It provides a consistent interface for accessing the number of probes across different classes.
#'
#' @param x An object from which to retrieve the number of probes.
#'
#' @return The number of probes as an integer.
#'
#' @export
setGeneric("nbprobes", function(x) standardGeneric("nbprobes"))

#' Method for retrieving the number of probes from a "Betas" object
#'
#' This method retrieves the number of probes from a "Betas" object.
#'
#' @param x A "Betas" object from which to retrieve the number of probes.
#'
#' @return The number of probes in the "Betas" object as an integer.
#'
#' @export
setMethod("nbprobes",signature("Betas"),definition=function(x) nrow(assays(x)$betas) )

#' Method for retrieving the number of probes from an "Mvals" object
#'
#' This method retrieves the number of probes from an "Mvals" object.
#'
#' @param x An "Mvals" object from which to retrieve the number of probes.
#'
#' @return The number of probes in the "Mvals" object as an integer.
#'
#' @export
setMethod("nbprobes",signature("Mvals"),definition=function(x) nrow(assays(x)$mvals) )


#' Create a new "Betas" object from beta values and a sample sheet
#'
#' This function creates a new "Betas" object from beta values and a sample sheet.
# It combines the provided data to create a "Betas" object with the appropriate structure.
#'
#' @param betas A matrix of beta values.
#' @param sampleSheet A data frame containing sample information.
#' @param na A numeric value, between 0 and 1, a probe is MASKED if its NA proportion is superior to na
#'
#' @return A "Betas" object containing the specified beta values, sample information, and MASK information as rowData.
#' 
#' @export
newBetas<-function(betas, sampleSheet, na) {
  
  # Check parameters
  assertthat::assert_that(is.matrix(betas), msg="betas must be a matrix of numeric values.")
  assertthat::assert_that(na >= 0 && na <= 1, msg="na must be between 0 and 1.")

  # check sampleSheet
  sampleSheet <- formatSampleSheet(sampleSheet) %>% S4Vectors::DataFrame()

  # Check for missing barcodes
  missing_in_betas <- setdiff(sampleSheet$barcode, colnames(betas))
  assertthat::assert_that( length(missing_in_betas)==0, msg=cat("Barcodes provided in sampleSheet are missing : ", missing_in_betas, sep=" ") )

  # If barcodes in betas are not in the sampleSheet, they are removed.
  missing_in_sampleSheet <- setdiff(colnames(betas), sampleSheet$barcode)
  if (length(missing_in_sampleSheet) > 0) {
    warning("The following barcodes are not in sampleSheet$barcode and will be removed from the analysis : ", paste(missing_in_sampleSheet, sep=","), "\n")
  }

  # Reorder sentrix
  rownames(sampleSheet)=sampleSheet$barcode
  samp_order<-dplyr::intersect(colnames(betas),sampleSheet$barcode)
  sampleSheet=sampleSheet[samp_order,]
  betas=betas[,samp_order]
  
  # create new Beta
  meth <- SummarizedExperiment( assays=list(betas=betas), colData=sampleSheet, rowData = getMask(betas) )
  meth <- as(meth,"Betas")
  metadata(meth)$plateform=getPlateform(betas)
  
  # create mask for NA
  probes <- CpGNAexcl(assays(meth)$betas,nalimit=na)
  rowData(meth)["Na"]=FALSE
  rowData(meth)["Na"][probes,]<-TRUE
  return(meth)
}

#' Compute M-values from a "Betas" object
#'
#' This method computes M-values from a "Betas" object, which represents beta values.
#'
#' @param x A "Betas" object from which to compute M-values.
#' @param grp A character vector specifying the group variable for differential analysis.
#' @param sva A character vector specifying the surrogate variable analysis (SVA) variable to include in the model. If NULL, SVA is skipped.
#' @param win A logical indicating whether Winsorization should be applied to the beta values. Default is TRUE.
#' @param sex A logical indicating whether sex chromosome data should be included. Default is FALSE.
#' @param sentrix2remove A list of Sentrix IDs to remove from the analysis. Default is an empty list.
#'
#' @return A "Mvals" object containing computed M-values.
#'
#' @export
setMethod("getMvals", signature("Betas"), definition = function(x, grp="grp", sva=NULL, win=TRUE, sex=FALSE, sentrix2remove=list() ) {  
  
  # Check parameters 
  assertthat::assert_that(class(x)=="Betas", msg="meth must be of class Betas")
  assertthat::assert_that(tolower(grp) %in% variables(x), msg="grp must be a column in the sampleSheet of betas")
  assertthat::assert_that( is.null(sva) || grepl("^~", sva), msg="sva is a model string an must start by ~")
  
  # get masked betas and remove probes with too many samples
  betas <- methylkey::getBetas(x, mask=FALSE, na=FALSE, sex=sex)
  
  # remove samples in sentrix2remove
  sel<-which(!colData(x)$barcode %in% sentrix2remove)
  pdata <- colData(x)[sel,]
  betas <- betas[,sel]
  
  # remove remaining na whith mean method
  group <- pdata %>% dplyr::pull(tolower(grp))
  betas <- replaceByMean( betas, groups = group )
  
  # winsorize betas
  if (win) {
    print("Winsorize data")
    betas <- DescTools::Winsorize(betas)
  }
  
  # calculate deltabetas
  print("Calculate deltabetas")
  deltaBetas=getDeltaBetas2(betas,group)
  
  # convert to mvalues
  print("Convert to mvalues")
  betas <- beta2m(betas)
  
  # batch correction with sva
  if (!is.null(sva)) {
    print("Batch correction")
    betas <- bc_sva(betas, pdata, sva)
  }
  
  mvals<-SummarizedExperiment( assays=list(mvals=betas), colData=pdata, rowData = deltaBetas )
  mvals <- as(mvals,"Mvals")
  metadata(mvals)=metadata(x)
  metadata(mvals)$grp=grp
  metadata(mvals)$sva=sva
  metadata(mvals)$winsorized=win
  metadata(mvals)$sex=sex
  
  return(mvals)
})

#' Determine the methylation array platform from the number of rows in a matrix
#'
#' This function attempts to determine the methylation array platform (e.g., IlluminaHumanMethylationEPIC, IlluminaHumanMethylation450k, IlluminaHumanMethylation27k, IlluminaMouseMethylation285k) based on the number of rows in a matrix.
# If the matrix size does not match any known platform, it returns "unknown."
#'
#' @param mat A matrix for which to determine the platform.
#'
#' @return The identified platform as a character string, or "unknown" if the platform cannot be determined.
#'
#' @export
getPlateform<-function(mat){
  
  p="unknown"
  if (is.null(mat)) return("unknown")
  if ( nrow(mat) == 937690){ return("IlluminaHumanMethylationEPICv2") }
  if ( nrow(mat) == 866553){ return("IlluminaHumanMethylationEPIC") }
  if ( nrow(mat) == 486427){ return("IlluminaHumanMethylation450k") }
  if ( nrow(mat) == 27578 ){ return("IlluminaHumanMethylation27k") }
  if ( nrow(mat) == 296070 ){ return("IlluminaMouseMethylation285k") }
  if ( nrow(mat) <= 937690){ p="IlluminaHumanMethylationEPICv2" }
  if ( nrow(mat) <= 866553){ p="IlluminaHumanMethylationEPIC" }
  if ( nrow(mat) <= 486427){ p="IlluminaHumanMethylation450k" }
  if ( nrow(mat) <= 296070){ p="IlluminaMouseMethylation285k" }
  if ( nrow(mat) <= 27578){ p="IlluminaHumanMethylation27k" }
  if ( nrow(mat) == 0){ p="unknown" }
  return(p)
}

#' Retrieve a masking matrix based on the Illumina methylation platform
#'
#' This function retrieves a masking matrix based on the Illumina methylation platform specified. The masking matrix is used for filtering and annotation of CpG probes.
#'
#' @param mat A data matrix representing methylation data.
#'
#' @return A masking matrix with relevant information for filtering CpG probes.
#'
#' @export
getMask<-function(mat){
  
  MASK=NULL
  plateform<-getPlateform(mat)
  print(plateform)
  if(plateform=="unknown") return(MASK)
  if (plateform=="IlluminaHumanMethylationEPICv2"){
    MASK<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.mask.tsv.gz", show_col_types = FALSE)
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz", show_col_types = FALSE)
  }
  if (plateform=="IlluminaHumanMethylationEPIC"){
    MASK<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/archive/202209/EPIC.hg38.mask.tsv.gz", show_col_types = FALSE) %>% 
      dplyr::rename(Probe_ID=probeID, M_general=MASK_general)
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz", show_col_types = FALSE)
  }
  if (plateform=="IlluminaHumanMethylationEPIC+"){
    MASK<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.mask.tsv.gz", show_col_types = FALSE)
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.manifest.tsv.gz", show_col_types = FALSE)
  }
  if (plateform=="IlluminaHumanMethylation450k"){
    MASK<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/archive/202209/HM450.hg38.mask.tsv.gz", show_col_types = FALSE) %>%
      dplyr::rename(Probe_ID=probeID, M_general=MASK_general)
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz", show_col_types = FALSE)
  }
  if(plateform=="IlluminaMouseMethylation285k"){
    MASK<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/MM285.mm10.mask.tsv.gz", show_col_types = FALSE)
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/MM285.mm10.manifest.tsv.gz", show_col_types = FALSE)
  }
  MASK=right_join(MASK,manifest[c("Probe_ID","CpG_chrm","CpG_beg","CpG_end")],by="Probe_ID") %>% tibble() %>%
    mutate(MASK_chrm=ifelse(CpG_chrm=="chrX" | CpG_chrm=="chrY",TRUE,FALSE))
  
  MASK$rs=ifelse(grepl("rs", MASK$Probe_ID) ,TRUE,FALSE)
  MASK$ctl=ifelse(grepl("ctl", MASK$Probe_ID),TRUE,FALSE)
  MASK<-MASK[ match(rownames(mat),MASK$Probe_ID,nomatch = NA), ]
  
  MASK$M_general <- ( MASK$M_general | MASK$rs | MASK$ctl )
  MASK$M_general[is.na(MASK$M_general)] <- FALSE
  
  return(MASK)
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
getManifest<-function(plateform){
  
  manifest=NULL
  
  if(plateform=="IlluminaHumanMethylation450k"){
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/HM450.hg38.manifest.tsv.gz")
  }
  
  if(plateform=="IlluminaHumanMethylationEPIC"){
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC/EPIC.hg38.manifest.tsv.gz")
  }
  
  if(plateform=="IlluminaHumanMethylationEPIC+"){
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPIC+/EPIC+.hg38.manifest.tsv.gz")
  }
  
  if(plateform=="IlluminaHumanMethylationEPICv2"){
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/EPICv2/EPICv2.hg38.manifest.tsv.gz")
  }
  
  if(plateform=="IlluminaMouseMethylation285k"){
    manifest<-readr::read_tsv("https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/MM285/MM285.mm10.manifest.tsv.gz")
  }
  
  manifest<-manifest %>% 
    dplyr::rename(chr=CpG_chrm, pos=CpG_beg, strand=mapYD_A) %>% 
    dplyr::select(Probe_ID,chr,pos,strand) %>%
    mutate(pos=pos+1) %>% mutate(strand = ifelse(is.na(strand), NA, ifelse(strand == "f", "+", "-")))
  
  return(manifest)
}