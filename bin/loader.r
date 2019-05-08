#!/usr/bin/Rscript
#####################################################################################
#
# Title  : loader.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################

loader<-function(){
  
  suppressPackageStartupMessages(library(wateRmelon))
  suppressPackageStartupMessages(library(data.table))
  
  #######################
  #1- Loading data
  print("loading data")
  pdata<-read.table(pdata, sep=separator, header=T)
  
  #convert variables to factor
  for (variable in colnames(pdata)){ 
  	if (length(levels(factor(pdata[,variable]))) <= nrow(pdata)/8) 
  	{ 
  		pdata[,variable]<-as.factor(pdata[,variable]) 
  	}  
  }
  
  #check samples
  tryCatch({
    samples<-pdata[samples]
  }, error = function(err) {
    stop(paste0("Check if '", samples, "' is a column in your pdata. You can specify the column's name for samples names with --samples "))
  })
  if(length(unique(samples))!=length(samples)){ 
    stop("Samples names must be uniques !")
  }
  
  #check barcode
  if (pipeline != "rnbeads"){
    tryCatch({
  	  barcode<-colnames(pdata[barcode])
    }, error = function(err) {
      stop(paste0("Check if '", barcode, "' is a column in your pdata. You can specify the column's name for barcodes with --barcode "))
    })
  }
  
  #check groups
  for (group in groups) {
    tryCatch({
      if ( !is.factor(pdata[,group] )){
        groups<-groups[-which(groups==group)]
        print(paste0("warning : ", group, " removed from groups because it is numeric." )) 
      }else{
        dir.create(paste(out, group, sep="/"))
      }
    }, error = function(err) {
      stop(paste0("Check if '", group, "' is a column in your pdata. You can specify the columns' names for groups with --groups "))
    })
  }
  groups<-colnames(pdata[groups])
  
  #default filters
  if (filters=="default"){
    filters<-paste0(datadir,c("Crossreactive_probes_EPIC.csv","SNP_EPIC.csv","SNP_EPIC_single_base.csv"))
  }
  
  
  #######################
  #2- Run pipeline
  if (pipeline=="minfi"){ 
  	source(paste0(path,"/pipeline_minfi.r"));
  	pdata[,barcode]<-as.character(pdata[,barcode]) 
  	pdata$Basename<-pdata[,barcode]
  }
  if (pipeline=="methylumi"){
   	source(paste0(path,"/pipeline_methylumi.r"))
  	if(normalize=="funnorm"){normalize="bmiq"} # default for methylumi is bmiq
  }
  if (pipeline=="rnbeads"){
   	source(paste0(path,"/pipeline_rnbeads.r"))
  	#pdata$barcode=as.character(samples)
  }
  
  #readmeth pipeline change depending the loaded pipeline
  data<-suppressPackageStartupMessages( readmeth(pdata=pdata, idat=idat, samples=unlist(samples), groups=groups, normalize=normalize, filters=filters, nalimit=nalimit, out=out, cell=cell) )
  
  regions=data$regions
  platform=data$platform
  pdata=data$pdata
  betas=data$betas
  nbprobes1=data$nbprobes1
  nbprobes2=data$nbprobes2
  filteredNAprobes=data$filteredNAprobes
  filteredFromList=data$filteredFromList
  
  #######################
  #3- PCA
  print("running PCA")
  save.image(file="debug.rdata")
  pca<-suppressPackageStartupMessages( makepca( betas, pdata, out, colnames(pdata), nPC=10 ) )
  pvalue<-pca$pvalue
  qvalue<-pca$qvalue
  jpeg(paste0(out, "/pca_contributions_0001.jpg"))
  plot(pca$contrib)
  dev.off()
  
  #######################
  #4- DeltaBetas
  #save.image(file=paste0(out, "/debug.rdata") )
  print("Calculate Delta Betas")
  deltab<-""
  deltab<-getDeltaBetas(betas,pdata[,groups[1]])
  
  #######################
  #5-violin plot
  if (violin){
  	print("Drawing violin plot")
  	save.image(file=paste0(out, "/debug.rdata") )
  	violin_plot(pdata=pdata, betas=betas, group=groups[1], samples=colnames(samples), platform=platform, path=path, out=out, genome=genome, regions=regions)
  }
  
  #######################
  #6- mvalues
  mval=beta2m(betas)
  mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
  # Remove a probe when sum of betas is 0.
  iszero<-which(rowSums(mval) == 0 )
  if(length(iszero) > 0){
  	mval<-mval[-which(rowSums(mval) == 0 ),]
  }
  
  #######################
  #6- Save and Create report
  print("Save results")
  write.table(pdata, file=paste0(out, "/pdata.txt"), row.names=F, sep="\t")
  write.table(betas, file=paste0(out, "/betas.txt"), row.names=T, sep="\t")
  write.table(deltab, file=paste0(out, "/deltabetas.txt"), row.names=T, sep="\t")
  write.table(pvalue, file=paste0(out, "/pvalue_0001.txt"), row.names=T, sep="\t")
  write.table(qvalue, file=paste0(out, "/qvalue_0001.txt"), row.names=T, sep="\t")
  
  save(betas, pdata, groups, samples, platform, genome, regions, mval, file=paste0(out, "/meth.rdata") )
  
  print("Create report")
  library(rmarkdown)
  if(pipeline == "rnbeads"){
    file.copy(paste0(path,"/rnbeads.Rmd"),"rnbeads.Rmd")
    render("rnbeads.Rmd", params = list(project=out, groups=groups ))
    file.remove("rnbeads.Rmd")
  }else{
    render(paste0(path,"/loader.Rmd"), params = list(project=out, groups=groups, workdir=getwd() ), knit_root_dir=getwd() )
  }
  
  
  
}




