#!/usr/bin/Rscript
#####################################################################################
#
# Title		: dmps.r
# Author	: CahaisV@iarc.fr
# Date		: 24/07/2018
# Last Update	: 24/07/2018
#
#####################################################################################


plot_dmrs<-function(){

  library(RColorBrewer)
  library(GenomicRanges)
  library(DMRcate)
  
  ########################
  if (platform=="sequencing"){ }
  if (platform=="IlluminaHumanMethylationEPIC"){ arraytype="EPIC"; require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) }
  if (platform=="IlluminaHumanMethylation450k"){ arraytype="450K"; require(IlluminaHumanMethylation450kanno.ilmn12.hg19) }
  
  #######################
  #3- colors
  lvl<-levels(as.factor(pdata[,groups[1]]))
  myColors <- brewer.pal(max(length(lvl),3),"Set1")
  names(myColors)<-lvl
  myColors<-myColors[as.character(pdata[,groups[1]])]
  
  #######################
  #4- dmrs ranges
  gdmr<-makeGRangesFromDataFrame(table)
  gdmr$dmrs<-table$dmrs
  
  ######################
  #dmr plot
  dmrplots<-list()
  for (dmr in dmrs){
  	if (dmr %in% table$dmrs){
  
  		dmrplots[dmr]<-mkPlotDMR(which(table$dmrs == dmr))
  	}
  }
  
  max<-min(length(table$dmrs),max)
  for (dmr in 1:max){
  	
    dmrname=table$dmrs[dmr]
    print(dmrname)
    img<-paste0(out, "/", model, "/dmrs/", dmrname, ".dmrplot.jpg")
    jpeg(img, width=800, height=800)
    DMR.plot(ranges=gdmr, dmr=dmr, CpGs=betas, phen.col=myColors, genome=genome, pch=20, toscale=T, plotmedians=T, arraytype=arraytype)
    dev.off()	
  }
  
  #######################
  #3-Create html report
  print("Create report")
  library(rmarkdown)
  render(paste0(path,"/plot_dmrs.Rmd"), params = list(project=out, model=model) )
  
}




