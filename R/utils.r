#' getDeltaBetas
#' 
#' calculate Delta betas between two groups
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

#' Calculate delta betas between all groups
#'
#' This function calculates delta betas (changes in methylation levels) between different groups based on methylation beta values. It computes the difference in methylation levels for all possible combinations of groups.
#'
#' @param betas A matrix of methylation beta values, where rows represent CpG sites and columns represent samples.
#' @param group A factor or character vector specifying group membership for each sample.
#'
#' @return A data frame containing delta betas for each pair of groups.
#'
#' @export
getDeltaBetas2 <- function(betas, group) {
  # Get unique levels from the 'group' factor
  levels_list <- levels(group)
  
  # Create an empty data frame to store the deltaBetas for each pair
  deltaBetas_df <- list()
  
  # Calculate deltaBetas for all combinations of levels
  for (i in 1:length(levels_list)) {
    for (j in 1:length(levels_list)) {
      case <- levels_list[i]
      control <- levels_list[j]
      if(case!=control){
        
        betas_<-betas[,group %in% c(case,control)]
        group_<-group[group %in% c(case,control)]
        deltaBetas <- rowMeans(betas_[,group==case]) - rowMeans(betas_[,group==control])
        
        deltaBetas_df[[paste0(case,"_vs_",control)]] <- deltaBetas * 100
      }
    }
  }
  
  return(data.frame(deltaBetas_df))
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
CpGexcl<-function(filters=NULL, plateform="IlluminaHumanMethylationEPIC", SNP=TRUE, CROSS=TRUE, XY=TRUE){
  
  probes=c()
  if(is.null(filters)){
    # from list
    probes=getDefaultProbeList(plateform, SNP, CROSS, XY)
  }
  else{
    # from files
    for (file in filters){
      if( !file.exists(file) & !RCurl::url.exists(file) ){ stop(paste0(file, " not found")) }
      probes <- unique( c(probes, readr::read_tsv(file)[[1]] ))
    }
  }
  
  return( probes )
}

#' makeSummarizedExperimentFromBetas
#'
#' create a summarized experiment
#'
#' @param betas betas array
#' @param sampleSheet data.frame
#'
#' @return SummarizedExperiment
#'
#' @export
#'
makeSummarizedExperimentFromBetas<-function(betas,sampleSheet){
  manifest=methylkey::getAnnotedManifest(methylkey::getPlateform(betas)) %>%
    filter(probeID %in% rownames(betas)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, start.field="pos", end.field="pos" )
  se<-SummarizedExperiment(assays=list(betas=betas[manifest$probeID,]), rowRanges=manifest, colData=sampleSheet)
  return(se)
}

