#####################################################################################
#
# Title  : tabular.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################




#######################################################
#        PIPELINE TO READ BETAS_1 FROM RNBEADS        #
#######################################################
#

readmeth<-function(pdata, idat , samples, groups, out, normalize="funnorm", cell="None"){
	
	print("pipeline for rnbeads data")

	#1-read idat files
	betas_1<- read.table(gzfile(idat), header=T, sep=",")
	regions<-betas_1[,1:5]
	betas<-as.matrix(betas_1[,6:ncol(betas_1)])
	colnames(betas)<-gsub("X","",colnames(betas))
	rownames(betas)<-regions$ID
	platform<-"RRBS"

	regions<-regions[,c(2,3,4,1)]
	regions<-cbind(regions,rep(".",nrow(regions)),rep("*",nrow(regions)) )
	write.table(regions, file=paste0(out,"/regions.bed"), row.names=F, col.names=F, quote=F, sep="\t" )
	regions<-paste0(out,"/regions.bed")

	return( list(betas=betas, pdata=pdata, platform=platform,regions=regions ) )

}
