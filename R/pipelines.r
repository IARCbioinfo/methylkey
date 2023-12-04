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
#' @import sesame
#' @import sesameData
#' @import wateRmelon
#'
#' @export
sesame2Betas<-function( idat=NULL, prep = "QCDPB", sampleSheet=NULL, na=0.2, ncore=4 ){
  
  require(sesame)
  require(sesameData)
  require(wateRmelon)
  require(purrr)
  
  #print("Warning loading data from EBV analysis")
  assertthat::assert_that( dir.exists(idat), msg=paste0("Directory not found : ", idat) )
  assertthat::assert_that( is.data.frame(sampleSheet), msg="Please provide a sampleSheet" )
  assertthat::assert_that(na >= 0 && na <= 1, msg="na must be between 0 and 1.")
  
  sdfs <- openSesame(idat, prep = prep, func = NULL, BPPARAM=BiocParallel::MulticoreParam(ncore))
  
  betas=openSesame(sdfs, prep = "", func = sesame::getBetas, mask = TRUE, BPPARAM=BiocParallel::MulticoreParam(ncore))
  
  sampleSheet<-formatSampleSheet(sampleSheet)
  inferedsex<-sapply(sdfs, function(sdf){ inferSex(sdf, platform = "EPIC", verbose = FALSE) }) %>%
    data.frame() %>% tibble::rownames_to_column("barcode") %>% dplyr::rename(inferedsex=2)
  sampleSheet<-merge(sampleSheet, inferedsex, by="barcode") %>% dplyr::relocate(samples)
  
  age<-wateRmelon::agep(betas, method='all') %>% tibble::rownames_to_column("barcode")
  sampleSheet<-merge(sampleSheet,age, by="barcode") %>% dplyr::relocate(samples)
  
  meth<-newBetas( betas, sampleSheet, na )
  metadata(meth)$betas.pipeline.name="sesame"
  
  qcs = openSesame(sdfs, prep="", func=sesameQC_calcStats, BPPARAM=BiocParallel::MulticoreParam(ncore))
  metadata(meth)$qcs = purrr::map( qcs, ~sesameQC_getStats(.) ) %>% bind_rows(.id="name")
  
  stx<-gsub("_.*","",names(sdfs))
  lapply(unique(stx), function(x){ plot_Channels(sdfs[stx==x]) })
  #metadata(meth)$RedGrnQQ=lapply(names(sdfs), function(x){ sesameQC_plotRedGrnQQ(sdfs[[x]]) })
  
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
#' @import minfi
#' @import wateRmelon
#'
#' @export
minfi2Betas<-function( idat=NULL, sampleSheet=NULL, na=0.2, compositeCellType="", pval=0.2 ){
  
  require(minfi)
  require(wateRmelon)
  message("001")
  compositeCellType_=c("Blood","CordBloodCombined","BloodExtended","CordBlood","CordBloodNorway","CordTissueAndBlood","DLPFC")
  cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu")
  if(compositeCellType=="BloodExtended"){
    cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv","CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")
  }
  message("002")
  assertthat::assert_that( dir.exists(idat), msg=paste0("Directory not found : ", idat) )
  assertthat::assert_that( "Basename" %in% colnames(sampleSheet), msg="Basename not in sampleSheet" )
  assertthat::assert_that( is.data.frame(sampleSheet), msg="Please provide a sampleSheet" )
  assertthat::assert_that(na >= 0 && na <= 1, msg="na must be between 0 and 1.")
  
  RGset<-read.metharray.exp(base = idat, targets=sampleSheet, force=TRUE)
  message("003")
  GMsetEx <- minfi::mapToGenome(RGset) 
  estSex  <- minfi::getSex(GMsetEx)
  sampleSheet$inferedsex<-estSex$predictedSex
  message("004")
  MSet <- minfi::preprocessRaw(RGset)
  MSet <- minfi::fixMethOutliers(MSet)
  qcs<-minfi::getQC(MSet)
  #Normalization
  betas <- minfi::getBeta(RGset)
  isna<-is.na(betas)
  MSet<-minfi::preprocessFunnorm(RGset, sex=estSex$predictedSex)
  betas <- minfi::getBeta(MSet)
  #After normalization NA values are replace by values close to 0. This restore the NA status.
  isna<-isna[ match(rownames(betas), rownames(isna)), ]
  betas[ which(isna) ]<-NA 
  message("005")
  pvalues <-minfi::detectionP(RGset)
  betas[ pvalues[rownames(betas), ] > pval ] <- NA
  message( paste0( "Low quality probes :", sum(pvalues > pval), " low quality probes replaced by NA"  ) )
  if (compositeCellType %in% compositeCellType_){
      #require(FlowSorted.Blood.EPIC)
    #require(FlowSorted.BloodExtended.EPIC)
    # FlowSorted.BloodExtended.EPIC <- libraryDataGet('FlowSorted.BloodExtended.EPIC')
    # FlowSorted.BloodExtended.EPIC
    # library(FlowSorted.BloodExtended.EPIC)
    # cc<-estimateCellCounts2(RGset, "BloodExtended", 
    #                         cellTypes = c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv","CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg"),
    #                         probeSelect = "IDOL",
    #                         referencePlatform = "IlluminaHumanMethylationEPIC",
    #                         CustomCpGs = IDOLOptimizedCpGsBloodExtended)
    
    if(getPlateform(betas) == "IlluminaHumanMethylationEPIC"){
      message("006")
      cc<-estimateCellCounts2(RGset, compositeCellType, 
                              cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu"),
                              referencePlatform = getPlateform(betas) )
      sampleSheet<-cbind(sampleSheet,cc$prop)
      sampleSheet$nlr=sampleSheet$neu/sampleSheet$nk
    }else{
      message("007")
      cc<-estimateCellCounts(RGset, compositeCellType, 
                              cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu"),
                              referencePlatform = getPlateform(betas) )
      sampleSheet<-cbind(sampleSheet,cc)
      sampleSheet$nlr=sampleSheet$Neu/sampleSheet$NK
    }
    message("008")
    breaks <- c(-Inf, 0.7, 1, 2, 3, Inf)
    labels <- c("bad", "average", "good", "average", "bad")
    sampleSheet$epimmune <- cut(sampleSheet$nlr, breaks = breaks, labels = labels)
  }
  message("009")
  sampleSheet<-cbind(sampleSheet, wateRmelon::agep(betas, method='all'))
  message("010")
  meth<-newBetas( betas, sampleSheet, na )
  message("011")
  metadata(meth)$betas.pipeline.name="minfi"
  metadata(meth)$qcs <- qcs
  message("012")
  plot_Channels2(RGset)
  message("013")
  return(meth)
}




