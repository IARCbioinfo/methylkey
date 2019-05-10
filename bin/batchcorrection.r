#!/usr/bin/Rscript
#####################################################################################
#
# Title  : batchcorrection.r
# Author : CahaisV@iarc.fr
# Date   : 04/07/2018
#
#####################################################################################


batchcorrection<-function(){
  
  suppressPackageStartupMessages(library(wateRmelon))

  ######################
  #1-check and reformat variables.
  print("check variables")
  variables<-suppressWarnings( autoformat(variables) )
  variables<-colnames(pdata[variables])
  for (v in variables) {
  	pdata[ ,v ]<-gsub(" ",".",pdata[ ,v ] )
  	pdata[ ,v ]<-gsub("-","_",pdata[ ,v ] )
  }
  
  ######################
  #2-create model
  print("create model")
  formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
  design<-model.matrix(formula1,data=pdata)
  
  ######################
  #3-mvalues
  print("check mvalues")
  mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
  
  # Remove a probe when sum of betas is 0.
  iszero<-which(rowSums(mval) == 0 )
  if(length(iszero) > 0){
  	mval<-mval[-which(rowSums(mval) == 0 ),]
  }
  
  ######################
  #4-Batch correction
  if(correction == "sva"){
  
  	print("running sva")
  
    if(nrow(design) > ncol(mval)){ stop("Missing samples in betas ! ") }
    if(nrow(design) < ncol(mval)){ stop("Missing samples in pdata ! ") }
    
  	require(sva)
  	
  	if (length(variables)==1){ 
  		sva<-sva(mval,design)
  	}
  	else{
  		formula0 <- as.formula(paste("~ 1+", paste(variables[-1],collapse="+")))
  		model0<-model.matrix(formula0, data=pdata)
  		sva<-sva(mval,design,model0)
  	}
  	#design<-cbind(design,sva$sv)
  	mval = t(residuals(lm(t(mval)~sva$sv)))
  
  	formula1<-paste(formula1, collapse="")
  	correction=paste0("on SVA corrected mvalues with model ", formula1)
  }
  
  if(correction == "ssva"){
  	
  	print("running SmartSVA")
  
  	require(SmartSVA)
  
  	#estimate number of surrogates
  	formula <- as.formula(paste("t(betas) ~ " ,paste(variables,collapse="+")))
  	Y.r <- t(resid(lm(formula, data=pdata)))
  	n.sv <- EstDimRMT(Y.r, FALSE)$dim
  	
  	if (length(variables)==1){ 
  		ssva<-smartsva.cpp(mval, design, n.sv=n.sv)
  	}
  	else{
  		formula0 <- as.formula(paste("~ 1+", paste(variables[-1],collapse="+")))
  		model0<-model.matrix(formula0, data=pdata)
  		ssva<-smartsva.cpp(mval, design, model0, n.sv=n.sv)
  	}
  	#design<-cbind(design,sva$sv)
  	mval = t(residuals(lm(t(mval)~sva$sv)))
  	
  	formula0<-paste(formula1, collapse="")
  	correction=paste0("on SSVA corrected mvalues with model ", formula1)
  }
  
  
  if(correction == "combat"){
  
  	print("running combat")	
  	require(sva)
  	mval<-ComBat(dat=mval, batch=pdata[,batch], mod=design, par.prior=TRUE, prior.plots=TRUE)
  	#if (!exists("sva")){sva=NULL}
  	correction=paste0("on combat corrected mvalues for batch ", colnames(pdata)[batch] )
  }
  
  
  #######################
  #6-PCA
  print("running PCA")
  
  #files<-list.files(out)
  #order<-length(grep("pca_contributions",files))
  #id=gsub("\\.","",as.character((order+1)/1000))
  
  pca<-suppressPackageStartupMessages(  makepca( mval, pdata, out, colnames(pdata), nPC=10 ) )
  qvalue<-pca$qval
  pvalue<-pca$pval
  
  jpeg(paste0(out, "/pca_contributions_0002.jpg"))
  plot(pca$contrib)
  dev.off()
  
  #######################
  #7- Save and Create report
  print("save results")
  write.table(pvalue, file=paste0(out, "/pvalue_0002.txt"), row.names=T, sep="\t")
  write.table(qvalue, file=paste0(out, "/qvalue_0002.txt"), row.names=T, sep="\t")
  write.table(mval, file=paste0(out, "/mval_batchcorrected.txt"), row.names=T, sep="\t")
  save(betas, mval, pdata, groups, samples, design, platform, genome, regions, file=paste0(out, "/meth.rdata") )

  print("Create report")
  library(rmarkdown)
  file.copy(paste0(path,"/batchcorrection.Rmd"),"batchcorrection.Rmd")
  render("batchcorrection.Rmd", params = list(project=out, model=formula1 ))
  file.remove("batchcorrection.Rmd")

}






