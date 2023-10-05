#' Methylkey differential analysis
#'
#' Create a differential methylation analysis report in Rmarkdown.
#'
#' @param se summarizedExperiment object with betas values in assays and illumina manifest as Features
#' @param model model to use for the analysis (string) eg: "~group+gender"
#' @param modelSVA model to use for the batch correction (string) if different from model.
#' @param case name of case level
#' @param control name of control level
#' @param sva if TRUE, sva batch correction will be applied to mvalues
#' @param method method to use in regression (ls or robust)
#' @param niter number of iteration with robust method
#' @param ncore number of core use for regression analysis
#' @param qval qvalue cutoff for selecting significant dmps
#' @param pcutoff qvalue cutoff for selecting dmps used for dmr analysis
#' @param sex if FALSE remove XY probes for analysis
#'
#' @return MethylResultSet
#'
#' @export
#' 
methyldiff2<-function(se=NULL, model=NULL,case=NULL,control=NULL,method="ls", niter=50, ncore=2, fdr=0.05, pcutoff=0.2, maxgap=1000, genome="hg38"){
  
  #model="~Group";modelSVA=NULL;case="D19";control="D0";sva=FALSE;pca=TRUE;method="ls";niter=50;ncore=2;qval=0.05;pcutoff=0.2
  
  manifest = getManifest( metadata(se)$plateform )
  pdata = data.frame(colData(se))
  mval = getMvals(se)
    
  #relevel according to case and control
  grp_g<-strsplit(model,"~|\\+")[[1]][2]
  pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), case)
  pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), control)
  
  #pull first group and sample names into vectors
  grp_g <- pdata[,grp_g] %>% as.factor()
  samples <- pdata$samples
  
  #regression analysis
  regression<-m_regression(mval,pdata, model, method=method, niter=niter, ncore=ncore )
  deltaBetas<-rowData(methM)[ rownames(regression$table), paste0(case,"_vs_",control) ]
  dmps<-data.frame(probeID=rownames(regression$table), regression$table, deltabetas=deltaBetas) %>%
    dplyr::select("probeID","P.Value","adj.P.Val","t","Coefficient","Stdev","deltabetas","goodness") %>%
    tidyr::separate(probeID, sep="_", into=c("probeID","x"),remove=TRUE, fill="right") %>%
    dplyr::mutate(status=ifelse(deltabetas>0,"hyper","hypo")) %>%
    dplyr::left_join(manifest,by="probeID") %>%
    tidyr::unite("probeID", c(probeID,x),sep="_", na.rm=TRUE, remove=TRUE) %>%
    dplyr::arrange(desc(adj.P.Val))
  
  dmrcate_table<-searchDMR_dmrcate(dmps,fdr=fdr,pcutoff=pcutoff,maxgap=maxgap,genome=genome) %>%
    dplyr::rename(fdr="HMFDR")
  #dmrff_table<-searchDMR_dmrff(dmps, mval[dmps$probeID,], maxgap=1000) %>%
    #dplyr::rename(no.cpgs="n", fdr="p.adjust")
  #combp_table<-searchDMR_combp(dmps, maxgap=maxgap) %>%
    #dplyr::rename(no.cpgs="nprobe")
  ipdmr_table<-searchDMR_ipdmr(dmps, maxgap=maxgap) %>%
    dplyr::rename(no.cpgs="nprobe")
  
  mrs<-MethylResultSet(se,dmps=dmps,dmrs=list(
    dmrcate=dmrcate_table,
    #dmrff=dmrff_table,
    #combp=combp_table,
    ipdmr=ipdmr_table))
  mrs@metadata=se@metadata
  mrs@metadata$model=model
  mrs@metadata$model=model
  mrs@metadata$model.case=case
  mrs@metadata$model.control=control
  mrs@metadata$model.method=method
  mrs@metadata$model.method.niter=niter
  mrs@metadata$model.method.ncore=ncore
  mrs@metadata$model.lambda=regression$lambda
  mrs@metadata$dmrs.fdr=fdr
  mrs@metadata$dmrs.pcutoff=pcutoff
  mrs@metadata$dmrs.maxgap=maxgap
  
  return(mrs)
} 






