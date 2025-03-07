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
getDeltaBetas<-function(betas,group){
  ( rowMeans(betas[,group==1, drop = FALSE]) - rowMeans(betas[,group==0, drop = FALSE]) ) *100
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
        
        #betas_<-betas[,group %in% c(case,control)]
        #group_<-group[group %in% c(case,control)]
        deltaBetas <- rowMeans(betas[,group==case]) - rowMeans(betas[,group==control])
        
        deltaBetas_df[[paste0(case,"_vs_",control)]] <- deltaBetas * 100
      }
    }
  }
  
  return(data.frame(deltaBetas_df, check.names = FALSE))
}

#' beta2m
#' 
#' calculate mvalues
#' 
#' @param betas array of betas values
#' 
#' @return matrix
#' 
#' @export
#' 
beta2m <- function (betas) {
  return(log2(betas/(1 - betas)))
}

#' beta2m
#' 
#' calculate mvalues
#' 
#' @param mvalues array of M values
#' 
#' @return matrix
#' 
#' @export
#' 
m2beta <- function(mvalues) {
  return(2^mvalues / (1 + 2^mvalues))
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




#' Generic function to extract column data from an object
#'
#' @param x An object RGset or sdfs list.
#' @param outfile file output name
#'
#' @export
setGeneric("toGeoSubmission", function(x,outfile) standardGeneric("toGeoSubmission") )


#' generate files to submit data to Geo database
#'
#' @param RGset an RGset object
#' @param outfile file output name
#'
#' @export
#'
setMethod("toGeoSubmission", signature(x="RGChannelSet",outfile="character"), definition = function(x,outfile="rawdata2Geo.tsv"){
  
  Mset<-minfi::preprocessRaw(RGset)
  meth<-minfi::getMeth(Mset)
  colnames(meth)<-paste( minfi::pData(RGset)$Basename, "Methylated signal")
  unmeth<-minfi::getUnmeth(Mset)
  colnames(unmeth)<-paste( minfi::pData(RGset)$Basename, "Unmethylated signal")
  pval<-minfi::detectionP(RGset)
  colnames(pval)<-paste( minfi::pData(RGset)$Basename, "Detection Pval")
  
  MSI<-cbind(unmeth,meth,pval)
  n<-ncol(pval)+1
  m<-ncol(pval)*2+1
  cols<-c(1,n,m)+rep(0:(n-2),each=3)
  MSI<-MSI[ , cols ]
  readr::write_tsv(as_tibble(MSI), quote=FALSE, file=outfile)
  
})


#' generate files to submit data to Geo database
#'
#' @param Betas a Betas object
#' @param outfile file output name
#'
#' @export
#'
setMethod("toGeoSubmission", signature(x="Betas",outfile="character"), definition = function(x,outfile="betas2Geo.tsv"){
  
  betas<-methylkey::getBetas(x)
  MP<-merge(betas,pval,by="row.names")
  rownames(MP)<-MP[,1]
  MP<-MP[,-1]
  cols<-c(1,n)+rep(0:(n-2),each=2)
  MP<-MP[ , cols ]
  write.table(MSI,sep="\t",row.names = T, quote=FALSE, file="../MatrixProcessed.tsv")
  
})
