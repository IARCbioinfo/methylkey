#!/usr/bin/Rscript
#####################################################################################
#
# Title  : dmrcate.r
# Author : CahaisV@iarc.fr
# Last update : 19/07/2018
#
#####################################################################################


dmrcate_wrapper<-function(){

  library(limma)
  library(DMRcate)
  library(data.table)
  
  print(model)
  
  if ( !file.exists(paste0(out,"/",model))) { dir.create(paste0(out,"/",model)) }
  if ( !file.exists(paste0(out,"/",model,"/dmrs"))) { dir.create(paste0(out,"/",model,"/dmrs")) }

  ########################
  #1- dmrcate 
  
  cmtx <- makeContrasts(colnames(design)[2], levels= colnames(design))
  colnames(design)[1] <- "(Intercept)"
  
  if (method=="robust"){ niter=200 }
  
  if (platform=="IlluminaHumanMethylationEPIC" | platform=="IlluminaHumanMethylation450k"){

  	arannot=NULL
  	if (platform=="IlluminaHumanMethylationEPIC"){ arraytype="EPIC"; arannot=c(array = "IlluminaHumanMethylationEPIC", annotation ="ilm10b2.hg19") }
  	if (platform=="IlluminaHumanMethylation450k"){ arraytype="450K"; arannot=c(array = "IlluminaHumanMethylation450k", annotation ="ilm10b2.hg19") }
  
  	myannotation<-cpg.annotate(object=mval, analysis.type=type, design=design, coef=2, datatype="array", arraytype=arraytype, what="M", fdr=fdr, annotation=arannot,  method=method, maxiter=niter, cont.matrix=cmtx)	
  	#dmrcoutput<-dmrcate(myannotation,lambda=1000,C=2,pcutoff=pcutoff)
  	dmrcoutput<-dmrcate(object=myannotation, C=2)
  	results.ranges <- extractRanges(dmrcoutput, genome=genome)
  }

  if (platform=="RRBS" | platform=="sequencing"){
    
  	final<-cbind(toptable$t, toptable$chr, toptable$start, toptable$logFC, toptable$adj.P.Val)
  	colnames(final)<-c('stat', 'chr', 'pos', 'diff', 'fdr')
  	myannotation<-cpg.annotate( as.data.frame(final), datatype="sequencing", analysis.type="differential", method=method, maxiter=niter)	
  	#suppressWarnings(suppressMessages(dmrcoutput<-dmrcate(myannotation,lambda=1000,C=50,pcutoff=0.05)))
  	myannotation$CHR<-paste0("chr",myannotation$CHR)
  	myannotation$CHR<-sub("chr23", "chrX", myannotation$CHR )
  	myannotation$CHR<-sub("chr24", "chrY", myannotation$CHR )
  	dmrcoutput<-dmrcate(myannotation,C=50)
  	#dmrcoutput$results$coord<-paste("chr",dmrcoutput$results$coord,sep="")	
  	#$results$coord<-sub("chr23", "chrX", dmrcoutput$results$coord )	
  	#dmrcoutput$results$coord<-sub("chr24", "chrY", dmrcoutput$results$coord )	
  	results.ranges <- extractRanges(dmrcoutput, genome = genome)
  }
  nbdmrs<-length(results.ranges)
  print( paste0("you get ", nbdmrs, " DMRS"))
  if (nbdmrs == 0) stop( "\r sorry !")
  
  #############################################
  #3-annotation
  print("Running annnotation")
  
  dmrs_annotated<-mk_annotate(results.ranges,genome=genome)
  
  #overlapping.probes
  locs <- GRanges(myannotation$CHR, IRanges(myannotation$pos, myannotation$pos, names=myannotation$ID))
  overlap <- findOverlaps(dmrs_annotated,locs)
  overlap <- data.table( cbind( query=queryHits(overlap), subject=subjectHits(overlap) ) )
  dmrs_annotated$overlapping.probes<-overlap[ , paste( names(locs)[subject], collapse="," ) , by=query ]$V1
  
  #dmrs
  dmrs<-data.frame(dmrs_annotated, row.names=NULL)
  colnames(dmrs)[1]<-"chr"
  dmrs$dmrs<-paste(dmrs$chr, dmrs$start, dmrs$end, sep=":")
  
  #remove duplicated lines
  dmrs<-dmrs[order(dmrs$annot.symbol),]
  dup<-duplicated(dmrs$dmrs)
  table<-dmrs[!dup,c(23,1:10,20,11,22)]
  
  #save files	
  dmrs<-dmrs[order(dmrs$minfdr),]
  table<-table[order(table$minfdr),]
  write.table(dmrs, file=paste(out, model, "dmrs", "annotated.txt",sep="/"), row.names=F, sep="\t")
  write.table(table, file=paste(out, model, "dmrs", "toptable.txt",sep="/"), row.names=F, sep="\t")
  
  ########################
  #barplots
  #save.image(file=paste0(out, "/debug.rdata"))
  
  p<-mk_barplot(annotated_regions=data.frame(dmrs_annotated), betafc=as.numeric(dmrs_annotated$meanbetafc), what="cpgi" )
  ggsave(file=paste(out, model, "dmrs", "cpgi_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300, units=)
  
  p<-mk_barplot(annotated_regions=data.frame(dmrs_annotated), betafc=as.numeric(dmrs_annotated$meanbetafc), what="genes" )
  ggsave(file=paste(out, model, "dmrs", "genes_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)
  
  ########################
  #11-Save Data
  print("Save data")
  save(betas, mval, pdata, groups, samples, design, platform, method, genome, regions, toptable, model, info, table, file=paste0(out, "/meth.rdata") )
  ########################
  #11-Create html report
  print("Create report")
  formula1<-paste(colnames(design)[-1], collapse=" + ")
  library(rmarkdown)
  file.copy(paste0(path,"/dmrs.Rmd"),"dmrs.Rmd")
  render("dmrs.Rmd", params = list(project=out, model=model, formula=formula1, method="dmrcate") )
  file.remove("dmrs.Rmd")
}










