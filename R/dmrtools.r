#' Search for Differentially Methylated Regions (DMRs) using DMRcate
#'
#' This function searches for Differentially Methylated Regions (DMRs) in DNA methylation data using the DMRcate package. It takes a set of CpG sites with associated statistical information and annotates them for DMR analysis.
#'
#' @param dmps A data frame containing CpG site information, including chromosome, position, strand, statistical test statistic, delta-betas, adjusted p-values, and probe IDs.
#' @param fdr The False Discovery Rate (FDR) threshold for identifying significant DMRs. Default is 0.05.
#' @param maxgap The maximum gap between neighboring CpGs to consider when defining DMRs. Default is 1000.
#' @param pcutoff The p-value threshold for filtering CpG sites before DMR analysis. Default is 0.2.
#' @param genome The genome build to use for genomic annotation. Default is "hg19."
#'
#' @return A data frame with information about identified DMRs, including chromosome, start and end positions, and associated CpG probe IDs.
#'
#'
#' @export
# searchDMR_dmrcate<-function(dmps, fdr=0.05, maxgap=1000,pcutoff=0.2,genome="hg19"){
#   require(DMRcate)
#   
#   annotated <- data.frame(chr=dmps$chr, start=dmps$pos, end=dmps$pos, strand=dmps$strand,
#                           stat=dmps$t, diff= dmps$deltabetas, ind.fdr=dmps$adj.P.Val, is.sig=(dmps$adj.P.Val<fdr) )
#   annotated<-GenomicRanges::makeGRangesFromDataFrame(annotated, keep.extra.columns=TRUE)
#   names(annotated)<-dmps$Probe_ID
#   myannotation <- new("CpGannotated", ranges=sort(annotated))
#   if( sum(is.na(myannotation@ranges$diff)) ){ myannotation@ranges$diff[ which(is.na(myannotation@ranges$diff)) ] <- 0 }
#   
#   dmrcoutput<- DMRcate::dmrcate(myannotation,C=2, pcutoff=0.2,lambda = maxgap)
#   table <- DMRcate::extractRanges(dmrcoutput, genome = genome)
#   
#   overlap <- GenomicRanges::findOverlaps(annotated,table,type="within")
#   final <- dmps[queryHits(overlap),]["Probe_ID"] %>%
#     bind_cols( as.data.frame(table)[subjectHits(overlap),] ) %>% 
#     dplyr::rename(chr="seqnames") %>%
#     dplyr::select(-strand) %>%
#     mutate(ID=paste0(chr, ":", start, "-", end) ) %>%
#     tidyr::nest(probes=Probe_ID)
#   
#   return(final)
# }


searchDMR_dmrcate<-function(dmps, fdr=0.05, maxgap=1000,pcutoff=0.2,genome="hg38",){
  
  if (!requireNamespace("DMRcate", quietly = TRUE)) {
    stop("Package 'DMRcate' is required for this function to work. Please install it.")
  }
  require(DMRcate)
  
  dmps <- dmps %>% dplyr::filter(!str_detect(chr, "_"))
  
  annotated <- data.frame(chr=dmps$chr, start=dmps$pos, end=dmps$pos, strand=dmps$strand,
                          stat=dmps$t, diff= dmps$deltabetas, ind.fdr=dmps$adj.P.Val, is.sig=(dmps$adj.P.Val<fdr) )
  annotated<-GenomicRanges::makeGRangesFromDataFrame(annotated, keep.extra.columns=TRUE)
  names(annotated)<-dmps$Probe_ID
  myannotation <- new("CpGannotated", ranges=sort(annotated))
  if( sum(is.na(myannotation@ranges$diff)) ){ myannotation@ranges$diff[ which(is.na(myannotation@ranges$diff)) ] <- 0 }
  
  dmrcoutput<- DMRcate::dmrcate(myannotation,C=2, pcutoff=pcutoff, lambda = maxgap)
  table <- DMRcate::extractRanges(dmrcoutput, genome = genome)
  overlap <- GenomicRanges::findOverlaps(annotated,table,type="within")
  #table <- as.data.frame(table)[subjectHits(overlap),c("seqnames","start","end","HMFDR","no.cpgs")]
  table <- as.data.frame(table)[subjectHits(overlap),c("seqnames","start","end","HMFDR","min_smoothed_fdr","no.cpgs")]
  #table$HMFDR <- table$min_smoothed_fdr

  dmps[queryHits(overlap),] %>% bind_cols( table )
}

#' Search for Differentially Methylated Regions (DMRs) using DMRff
#'
#' This function searches for Differentially Methylated Regions (DMRs) in DNA methylation data using the DMRff package. It takes a set of CpG sites with associated statistical information and methylation data and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information, including chromosome, position, statistical coefficient, standard deviation, p-value, and probe IDs.
#' @param betas A data frame containing methylation data, typically beta values, for the CpG sites.
#' @param maxgap The maximum gap between neighboring CpGs to consider when defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs, including chromosome, start and end positions, associated CpG probe IDs, and DMR tool name.
#'
#' @export
searchDMR_dmrff<-function(dmps, betas, maxgap=1000){
  
  if (!requireNamespace("dmrff", quietly = TRUE)) {
    stop("Package 'dmrff' is required for this function to work. Please install it.")
  }
  require(dmrff)
  
  dmrs <- dmrff(estimate=dmps$Coefficient,
                se=dmps$Stdev,
                p.value=dmps$P.Value,
                methylation=betas,
                chr=as.vector(dmps$chr),
                pos=dmps$pos,
                maxgap=maxgap,
                verbose=T)
  
  #some ranges are inverted
  dmrs<-dmrs %>% mutate(start = if_else(start > end, end, start), end = if_else(start > end, start, end))
  
  #find overlaps
  gr1<-GenomicRanges::makeGRangesFromDataFrame(dmps, start.field = "pos", end.field = "pos")
  gr2<-GenomicRanges::makeGRangesFromDataFrame(dmrs)
  overlap <- GenomicRanges::findOverlaps(gr1,gr2,type="within")
  
  dmrs<-dmrs %>% dplyr::select(-chr) # duplicated row with dmps
  final <- dmps[queryHits(overlap),] %>%
  bind_cols( dmrs[subjectHits(overlap),] ) %>%
  mutate(ID=paste0(chr, ":", start, "-", end), dmrtool="dmrff") %>%
  tidyr::nest(probes=Probe_ID)
  
  return(final)
}

#' Search for Differentially Methylated Regions (DMRs) using combp
#'
#' This function searches for Differentially Methylated Regions (DMRs) in DNA methylation data using the combp package. It takes a set of CpG sites with associated statistical information and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information, including chromosome, position, p-value, and probe IDs.
#' @param maxgap The maximum gap between neighboring CpGs to consider when defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs, including chromosome, start and end positions, associated CpG probe IDs, and DMR tool name.
#'
#' @export
searchDMR_combp<-function(dmps, maxgap=1000){
  
  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required for this function to work. Please install it.")
  }
  require(ENmix)
  
  data=data.frame(probe=rownames(dmps),p=dmps$P.Value,chr=dmps$chr,start=dmps$pos,end=dmps$pos)
  
  combp(data,dist.cutoff=1000,bin.size=310,seed=0.05, region_plot=FALSE,mht_plot=FALSE,nCores=10,verbose=TRUE)
  dmrs=readr::read_csv("resu_combp.csv")
  
  gr1<-GenomicRanges::makeGRangesFromDataFrame(dmps, start.field = "pos", end.field = "pos")
  gr2<-GenomicRanges::makeGRangesFromDataFrame(dmrs)
  overlap <- GenomicRanges::findOverlaps(gr1,gr2,type="within")
  
  dmrs<-dmrs %>% dplyr::select(-chr) # duplicated row with dmps
  final <- dmps[queryHits(overlap),] %>%
    bind_cols( dmrs[subjectHits(overlap),] ) %>%
    mutate(ID=paste0(chr, ":", start, "-", end), dmrtool="combp") %>%
    tidyr::nest(probes=Probe_ID)
  
  return(final)
}

#' Search for Differentially Methylated Regions (DMRs) using ipdmr
#'
#' This function searches for Differentially Methylated Regions (DMRs) in DNA methylation data using the ipdmr package. It takes a set of CpG sites with associated statistical information and performs DMR analysis.
#'
#' @param dmps A data frame containing CpG site information, including chromosome, position, p-value, and probe IDs.
#' @param maxgap The maximum gap between neighboring CpGs to consider when defining DMRs. Default is 1000.
#'
#' @return A data frame with information about identified DMRs, including chromosome, start and end positions, associated CpG probe IDs, and DMR tool name.
#'
#' @export
# searchDMR_ipdmr<-function(dmps, maxgap=1000){
#   require(ENmix)
#   
#   data=data.frame(probe=rownames(dmps),p=dmps$P.Value,chr=dmps$chr,start=dmps$pos,end=dmps$pos)
#   
#   ipdmr(data,dist.cutoff=1000,bin.size=310,seed=0.05, region_plot=FALSE,mht_plot=FALSE,verbose=TRUE)
#   dmrs=readr::read_csv("resu_ipdmr.csv")
#   
#   gr1<-GenomicRanges::makeGRangesFromDataFrame(dmps, start.field = "pos", end.field = "pos")
#   gr2<-GenomicRanges::makeGRangesFromDataFrame(dmrs)
#   seqlevelsStyle(gr2) <- "UCSC"
#   overlap <- GenomicRanges::findOverlaps(gr1,gr2,type="within")
#   
#   dmrs<-dmrs %>% dplyr::select(-chr) # duplicated row with dmps
#   final <- dmps[queryHits(overlap),] %>%
#     bind_cols( dmrs[subjectHits(overlap),] ) %>%
#     mutate(ID=paste0(chr, ":", start, "-", end), dmrtool="ipdmr") %>%
#     tidyr::nest(probes=Probe_ID)
#   
#   return(final)
# }

searchDMR_ipdmr<-function(dmps, maxgap=1000, bin.size=310, seed=0.05){
  
  if (!requireNamespace("ENmix", quietly = TRUE)) {
    stop("Package 'ENmix' is required for this function to work. Please install it.")
  }
  
  require(ENmix)
  
  data=data.frame(probe=dmps$Probe_ID,p=dmps$P.Value,chr=dmps$chr,start=dmps$pos,end=dmps$pos)
  data$p[is.na(data$p)]<-0.99999

  ipdmr<-data.frame(chr=character(),start=numeric(),end=numeric(),p=numeric(),fdr=numeric(),nprobe=numeric(),probe=character())
  write_csv(ipdmr,"resu_ipdmr.csv") # avoid to reload previous results if ipdmr find 0 dmrs.
  ipdmr(data, dist.cutoff=maxgap, bin.size=bin.size, seed=seed, region_plot=FALSE, mht_plot=FALSE, verbose=FALSE)
  dmrs=readr::read_csv("resu_ipdmr.csv")
  
  return(dmrs)
}
