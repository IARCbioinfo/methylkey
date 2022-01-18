#' getDeltaBetas
#' 
#' calculate Delta betas
#' 
#' @param betas array of betas values
#' @param group group
#' @param case case
#' @param control control
#' 
#' @return vector
#' 
#' @export
#' 
getDeltaBetas<-function(betas,group, case="TT", control="NT"){
  betas<-betas[,group %in% c(case,control)]
  group<-group[group %in% c(case,control)]
  deltaBetas <- rowMeans(betas[,group==case]) - rowMeans(betas[,group==control])
  return(deltaBetas)
}

#' beta2m
#' 
#' calculate mvalues
#' 
#' @param betas array of betas values
#' 
#' @return vector
#' 
#' @export
#' 
beta2m <- function (betas) {
  return(log2(betas/(1 - betas)))
}

#' getDefaultProbeList
#' 
#' get default list of probes to filter according to the plateform
#' 
#' @param filters file containing probe list
#' 
#' @return vector
#' 
getDefaultProbeList<-function(plateform="IlluminaHumanMethylationEPIC", SNP=TRUE, CROSS=TRUE, XY=FALSE){
  
  filters<-c()
  if (plateform=="IlluminaHumanMethylationEPIC"){
    if (CROSS){ filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/Crossreactive_probes_EPIC.csv" ) }
    if (SNP)  { filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/SNP_EPIC.csv" ) }
    if (SNP)  { filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/SNP_EPIC_single_base.csv" ) }
    if (XY)   { filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/Sex_EPIC.csv" ) }
  }
  
  if (plateform=="IlluminaHumanMethylation450k"){
    if (CROSS){ filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/Crossreactive_probes_450k_Chen.csv") }
    if (CROSS){ filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/Crossreactive_probes_450k_EGE.csv") }
    if (SNP)  { filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/SNP_All_races_5percent.csv") }
    if (XY)   { filters<- c(filters, "http://git.iarc.lan/EGE/methylkey/raw/master/data/Sex_Chr_SNPs.csv") }
  }
  
  return(filters)
}

#' CpGexcl
#' 
#' Filter probes from list of probes
#' 
#' @param filters file containing probe list
#' @param plateform plateform
#' 
#' @return vector
#' 
#' @export
#' 
CpGexcl<-function(filters=NULL, plateform="IlluminaHumanMethylationEPIC"){
  
  # from list
  if(is.null(filters)){
    filters=getDefaultProbeList(plateform)
  }
  
  probes=c()
  for (file in filters){
    if( !file.exists(file) & !RCurl::url.exists(file) ){ stop(paste0(file, " not found")) }
    probes <- unique( c(probes, read_tsv(file)[[1]] ))
  }
  
  return( probes )
}

#' getPlateform
#'
#' guess illumina plateform based on probes names and probes number
#'
#' @param betas betas array
#'
#' @return plateform
#'
#' @export
#'
getPlateform<-function(betas){
  
  if ( grepl("_",rownames(betas)[1]) ){ return("IlluminaMouseMethylation285k") }
  else if ( nrow(betas) > 500000){ return("IlluminaHumanMethylationEPIC") }
  else if ( nrow(betas) > 200000){ return("IlluminaHumanMethylation450k") }
  else if ( nrow(betas) > 10000 ){ return("IlluminaHumanMethylation27k") }
  else { stop("can't identify plateform !") }
}
