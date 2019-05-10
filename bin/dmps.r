#!/usr/bin/Rscript
#####################################################################################
#
# Title		: dmps.r
# Author	: CahaisV@iarc.fr, RenardC@iarc.fr
# Date		: 07/10/2014
# Last Update	: 05/07/2018
#
#####################################################################################


dmps<-function(){
  
  library(GetoptLong)
  library(limma)
  suppressPackageStartupMessages(library(qqman))
  
  #############################################
  #2-Prepare variables
  variables<-suppressWarnings( autoformat(variables) )
  variables<-colnames(pdata[variables])
  pdata[ ,variables[1] ]<-as.factor(pdata[ ,variables[1] ])
  for (v in variables) {
  	if (is.factor(v)){
  		pdata[ ,variables[v] ]<-make.names(pdata[ ,variables[v] ])
  	}
  }
  
  #############################################
  #3-Order the levels : control become intercept and case is compared to control.
  case<-make.names(case)
  control<-make.names(control)
  pdata[,variables[1]] <- relevel(pdata[,variables[1]], case)
  pdata[,variables[1]] <- relevel(pdata[,variables[1]], control)
  
  #############################################
  #4-create model
  formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
  if (formula!="None"){ formula1 <- as.formula(paste("~ " ,formula, collapse="") ) }
  model=paste0(case,"_vs_",control)
  if ( !file.exists(paste0(out,"/",model))) { dir.create(paste0(out,"/",model)) }
  if ( !file.exists(paste0(out,"/",model,"/dmps"))) { dir.create(paste0(out,"/",model,"/dmps")) }
  print(model)
  
  design<-model.matrix(formula1,data=pdata)
  colnames(design)<-make.names(colnames(design))
  cmtx <- makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )
  #############################################
  #5-remove cpgs when sum of its betas is 0
  #it can produce a BUG when it is the case
  iszero<-which(rowSums(mval) == 0 )
  if(length(iszero) > 0){
  	mval<-mval[-which(rowSums(mval) == 0 ),]
  }
  
  #############################################
  #6-run model
  print("running lmfit")
  fit<-lmFit(mval,design, pdata, ndups=1, method=method)
  rownames(cmtx)<-colnames(fit)
  fitContrasts=contrasts.fit(fit,cmtx)
  eb=eBayes(fitContrasts)
  lmPvals = eb$p.value
  chisq <- qchisq(1-eb$p.value,1)
  lambda<-median(chisq)/qchisq(0.5,1)
  #qqplot
  jpeg(paste(out, model, "dmps/qqplot.jpg", sep="/"))
  suppressWarnings( qq(lmPvals,main=paste("QQ plot: ","lambda=",lambda,sep="")) )
  dev.off()
  
  #############################################
  #7-get top dmps
  toptable<-topTable(eb, adjust="BH", number=Inf, p=fdr, sort.by="P")
  nbdmps=nrow(toptable)
  print( paste0("you get ", nbdmps, " DMPS"))
  if (nbdmps == 0) {
  	#save(betas, mval, ssva, pdata, groups, samples, design, platform, html, out, file=meth )
  	stop( "\r sorry !")
  }
  
  #############################################
  #8-delta betas
  print("Calculating delta betas")
  deltabetas<-getDeltaBetas(betas[rownames(toptable),],pdata[,variables[1]])
  toptable$deltaBetas<-deltabetas[,1]
  
  #############################################
  #9-annotation
  print("Running annnotation")
  
  #get regions
  if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
  if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
  if(regions == "") { stop( "\r Oups ! no region to annotate.") }
  
  #complete annotation
  cpg_annotated<-mk_annotate(regions,genome=genome)
  dm_annotated<-cpg_annotated[ cpg_annotated$name %in% rownames(toptable)]
  dm<-merge(toptable, dm_annotated, by.x="row.names", by.y="name")
  colnames(dm)[1]<-"cpg"
  colnames(dm)[9]<-"chr"
  dm$DM_Status[ as.numeric(dm$logFC) > 0 ] = "hyper"
  dm$DM_Status[ as.numeric(dm$logFC) < 0 ] = "hypo"
  write.table(dm, file=paste(out, model, "dmps", "annotated.txt",sep="/"), row.names=F, sep="\t")
  
  #toptable
  dm_annotated<-makeGRangesFromDataFrame(dm, keep.extra.columns=TRUE)
  dm<-dm[order(dm$annot.symbol),]
  toptable<-dm[!duplicated(dm$cpg),c(9:11,1,13,2:8,23)]
  toptable<-toptable[order(toptable$P.Value),]
  write.table(toptable, file=paste(out, model, "dmps", "toptable.txt",sep="/"), row.names=F, sep="\t")
  
  #############################################
  #10- Plots
  
  #circusplot
  #save.image(file="debug.rdata")
  
  print("circus plot")
  gbm<-makeGRangesFromDataFrame(toptable)
  gbm$meth.diff<-toptable$deltaBetas*100
  #suppressMessages ( p<-circusplot(gbm, genome=genome) )
  jpeg(paste(out, model, "dmps", "circleplot.jpg", sep="/") , width=800, height=800)
  #print(p)
  dev.off()
  
  #heatmap
  save.image("debug.RData")
  print("heatmap")
  require(NMF)
  if (nrow(toptable)<hsize){ hsize=nrow(toptable) }
  
  if (hsize>1) {
  	cpgs<-toptable$cpg[1:hsize]
  	symbols<-toptable$annot.symbol[1:hsize]
  	labrow<-paste0( cpgs, "(", symbols, ")")
  	jpeg(paste(out, model, "dmps", "heatmap.jpg",sep="/"), width=1600, height=1600)
  	aheatmap(betas[cpgs,], annCol=pdata[,groups], labRow=labrow)
  	dev.off()
  }
  
  ########################
  #11- barplots
  
  p<-mk_barplot(annotated_regions=dm_annotated, betafc=as.numeric(dm$logFC), what="cpgi" )
  ggsave(file=paste(out, model, "dmps", "cpgi_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300, units=)
  
  p<-mk_barplot(annotated_regions=dm_annotated, betafc=as.numeric(dm$logFC), what="genes" )
  ggsave(file=paste(out, model, "dmps", "genes_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)
  
  
  ########################
  #12- histone mark
  p<-mk_hmplot(dm, out, model, celllines)
  ggsave(file=paste(out, model, "dmps", "chromHMM_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)
  #save.image(file=paste0(out,"/debug.rdata"))
  
  ########################
  #MAPINFO
  print("create info file")
  
  info<-unique(data.frame(cpg_annotated)[,c(1,2,6)])
  info<-merge(info, eb$p.value, by.x="name", by.y="row.names")
  #info<-info[,c(1,2,3,8)]
  colnames(info)<-c("TargetID", "CHR", "MAPINFO", "Pval")
  #rownames(info)<-seq(1:nrow(info))
  
  formula1<-paste(formula1, collapse=" ")
  
  print("Save data")
  save(betas, mval, pdata, groups, samples, design, platform, method, genome, regions, toptable, model, info, file=paste0(out, "/meth.rdata") )
  
  print("Create report")
  library(rmarkdown)
  file.copy(paste0(path,"/dmps.Rmd"),"dmps.Rmd")
  render("dmps.Rmd", params = list(project=out, model=model, formula=formula1) )
  file.remove("dmps.Rmd")
}


