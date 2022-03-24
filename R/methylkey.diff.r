#' Methylkey differential analysis
#'
#' Create a differential methylation analysis report in Rmarkdown.
#'
#' @param model model to use for the analysis (string) eg: "~group+gender"
#' @param modelSVA model to use for the batch correction (string) if different from model.
#' @param case name of case level
#' @param control name of control level
#' @param betas matrix of beta values
#' @param pdata SampleSheet (dataframe)
#' @param sva if TRUE, sva batch correction will be applied to mvalues
#' @param method method to use in regression (ls or robust)
#' @param niter number of iteration with robust method
#' @param ncore number of core use for regression analysis
#' @param qval qvalue cutoff for selecting significant dmps
#' @param pcutoff qvalue cutoff for selecting dmps used for dmr analysis
#' @param output output folder
#' @param manifest annotation for th array
#' @param display number of dmps displayed in heatmap
#' @param plateform plateform type (IlluminaHumanMethylation450k, IlluminaHumanMethylationEPIC, IlluminaMouseMethylation285k)
#' @param genome genome assembly
#'
#' @return None
#'
#' @export
#' 
methyldiff<-function(model=NULL,modelSVA=NULL,case=NULL,control=NULL,betas=NULL,pdata=NULL,sva=TRUE,method="ls", niter=50,ncore=2,qval=0.05,pcutoff=0.2, output="methylkey", manifest=NULL, display=500,plateform=NULL,genome="hg19",level="###" ){
  
  #relevel according to case and control
  grp_g<-strsplit(model,"~|\\+")[[1]][2]
  pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), case)
  pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), control)
  
  #pull first group and sample names into vectors
  grp_g<- pdata %>% pull(grp_g) %>% as.factor()
  samples <- pdata %>% pull(samples)
  
  #identify plateform
  if(is.null(plateform)){
    plateform=getPlateform(betas)
    message(paste0("plateform=",plateform))
    print(plateform)
  }
  
  #annotated manifest
  if(is.null(manifest)){
    manifest <- getAnnotedManifest( plateform ) %>% 
      plyr::rename(replace = c(Name="probeID"), warn_missing = FALSE)
  }
  
  #correct default genome
  if (plateform=="IlluminaMouseMethylation285k") { genome=="mm10" }
  
  #calculate mvalues
  mval<-beta2m(betas)
  mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
  #iszero<-which(rowSums(mval) == 0 )
  
  cat(paste0("\n\n",level," Summary\n\n"))
  
  #logs
  print(paste0("model=",model))
  print(paste0("case=",case))
  print(paste0("control=",control))
  print(paste0("sva=",sva))
  print(paste0("method=",method))
  print(paste0("qval=",qval ,"  (Filter for DMPs and DMRs significance)"))
  print(paste0("pcutoff=",pcutoff, "  (Filter for selecting DMPs to use in DMR analysis)"))
  print(paste0("output=",output))
  
  #batch correction
  if(sva){
    cat(paste0("\n\n",level," Batch correction \n\n"))
    if(is.null(modelSVA)){ modelSVA=model }
    mval  <- bc_sva(mval,pdata,modelSVA)
  } else {
    cat(paste0("\n\n",level," PCA \n\n"))
  }

  #pca
  pca<-makepca( mval, pdata, nPC=9 )
  rownames(pca$pca$x)<-pdata$samples
  fig_sva<-factoextra::fviz_pca_ind(pca$pca,habillage=grp_g, addEllipses=TRUE, ellipse.level=0.95)
  print(fig_sva)
  print(pca$contrib)
  pca$pvalue %>% as.data.frame() %>% DT::datatable() %>% htmltools::tagList() %>% print()
  
  cat(paste0("\n\n",level," DMPs {.tabset}\n\n"))
  
  #regression analysis
  regression<-m_regression(mval,pdata, model, method=method, niter=niter, ncore=ncore )
  deltaBetas<-getDeltaBetas(betas[rownames(regression$table),],grp_g,case,control)*100
  table<-data.frame(probeID=rownames(regression$table), regression$table, deltabetas=deltaBetas) %>%
    dplyr::select("probeID","P.Value","adj.P.Val","t","Coefficient","Stdev","deltabetas") %>%
    separate(probeID, sep="_", into=c("probeID","x"),remove=TRUE, fill="right") %>%
    dplyr::mutate(status=ifelse(deltabetas>0,"hyper","hypo")) %>%
    dplyr::mutate(is.sig=adj.P.Val<qval) %>%
    merge(manifest, by="probeID") %>%
    tidyr::unite("probeID", c(probeID,x),sep="_", na.rm=TRUE, remove=TRUE) %>%
    dplyr::arrange(desc(adj.P.Val))
    
  write_tsv(table, file=paste0(output, "_", Sys.Date(), ".dmps.csv"))

  #qqplot
  cat(paste0("nb dmps=",regression$table %>% filter(adj.P.Val < qval) %>% nrow()))
  qqman::qq(regression$pvals,main=paste("QQ plot: ","lambda=", regression$lambda,sep=""))
  
  #toptable
  cat(paste0("\n\n#",level," TopTable (500 max) \n\n"))
  table %>% dplyr::filter(adj.P.Val < qval) %>% head(500) %>% DT::datatable() %>% htmltools::tagList() %>% print()
  
  #plots
  cat(paste0("\n\n#",level," DeltaBetas distribution\n\n"))
  deltaBetas<-deltaBetas[ regression$table$adj.P.Val< qval ]
  if (length(deltaBetas)>0) { hist(deltaBetas, col="chocolate" ) }

  cat(paste0("\n\n#",level," Volcano plot\n\n"))
  print(methylkey::volcano(table))
  
  cat(paste0("\n\n#",level," Manhattan plot\n\n"))
  print(methylkey::manhattan(table))
  
  cat(paste0("\n\n#",level," Barplots\n\n"))
  barplots_(plateform=plateform,table,manifest)
  
  #pathway analysis : enrichR
  if( (table %>% dplyr::filter(adj.P.Val < qval) %>% nrow()) > 3 ){
    
    cat(paste0("\n\n#",level," Heatmap\n\n"))
    NMF::aheatmap( betas[ table$probeID[1:min(100,nrow(table))] ,] , labCol=samples , annCol = grp_g, Rowv= FALSE )
  
    cat(paste0("\n\n#",level," EnrichR\n\n"))
    genes<-table %>% filter(adj.P.Val < qval) %>% dplyr::select(UCSC_RefGene_Name) %>% head(1000) %>% pull() %>% strsplit(";") %>% unlist() %>% unique() 
    enriched <- enrichR::enrichr(genes, dbs)
    for( db in dbs)
    {
      print(db)
      try(
        print(plotEnrich(enriched[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = db))
      )
    }
    
  }
  
  #DMRcate
  cat(paste0("\n\n",level," DMRs {.tabset}\n\n"))
  table     <- GenomicRanges::makeGRangesFromDataFrame(table, start.field = "pos", end.field= "pos", keep.extra.columns=TRUE)
  annotated <- data.frame(chr=GenomicRanges::seqnames(table), start=start(table), end=end(table), strand=strand(table),
                          stat=table$t, diff= table$deltabetas, ind.fdr=table$adj.P.Val, is.sig=(table$adj.P.Val<qval) )
  annotated<-GenomicRanges::makeGRangesFromDataFrame(annotated, keep.extra.columns=TRUE)
  names(annotated)<-table$probeID
  myannotation <- new("CpGannotated", ranges=sort(annotated))

  message("debug03")
  
  # fix issue when there is missing values in deltabetas.
  if( sum(is.na(myannotation@ranges$diff)) ){myannotation@ranges$diff[ which(is.na(myannotation@ranges$diff)) ] <- 0 }

  #for debug
  #myannotation2 <- cpg.annotate("array", mval, arraytype = "EPIC",analysis.type="differential", design=model.matrix(as.formula(model),data=pdata), coef=2, what = "M")
  dmrcoutput=NULL
  tryCatch(
    dmrcoutput<- DMRcate::dmrcate(myannotation,C=2, pcutoff=pcutoff)
  , error=function(e) { message(e) } )
  if(is.null(dmrcoutput)) return()
  
  table <- DMRcate::extractRanges(dmrcoutput, genome = genome)

  overlap <- GenomicRanges::findOverlaps(table,annotated)
  overlap <- cbind( query=queryHits(overlap), subject=names(annotated)[subjectHits(overlap)] )
  overlap <- overlap %>% data.frame() %>%
                mutate(query=as.numeric(query)) %>%
                group_by(query) %>% 
                summarise( overlapping.probes=paste0(subject,collapse = ",") )
  
  table <- as.data.frame(table)
  table$overlapping.sites<-overlap$overlapping.probes
  
  cat(paste0("\n\n#",level," TopTable (500 max)\n\n"))
  cat(paste0("nb dmrs=", table %>% filter(HMFDR < qval) %>% nrow()))
  table %>% filter(HMFDR < qval) %>% head(500) %>% DT::datatable() %>% htmltools::tagList() %>% print()
  
  write_tsv(table, file=paste0(output, "_", Sys.Date(), ".dmrs.csv"))
  
  #plots
  cat(paste0("\n\n#",level," DeltaBetas distribution\n\n"))
  meandiff<-table %>% filter( HMFDR< qval ) %>% pull(meandiff)
  if (length(meandiff)>0) { hist(meandiff, col="chocolate" ) }
  maxdiff<-table %>% filter( HMFDR< qval ) %>% pull(maxdiff)
  if (length(maxdiff)>0) { hist(maxdiff, col="chocolate" ) }
  
  cat(paste0("\n\n#",level," Volcano plot\n\n"))
  print(methylkey::volcano(table))
  
  cat(paste0("\n\n#",level," Manhattan plot\n\n"))
  print(methylkey::manhattan(table))
  
}
