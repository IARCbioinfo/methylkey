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

readmeth<-function(pdata, idat , samples, groups, normalize="funnorm", filters=filters, nalimit=nalimit, out, cell="None"){

	print("pipeline for rnbeads data")

	#1-read idat files
	betas_1<- read.table(gzfile(idat), header=T, sep=",")
	regions<-betas_1[,1:5]
	betas<-as.matrix(betas_1[,6:ncol(betas_1)])
	rownames(betas)<-regions$ID
	platform<-"RRBS"

	regions<-regions[,c(2,3,4,1)]
	regions<-cbind(regions,rep(".",nrow(regions)),rep("*",nrow(regions)) )
	write.table(regions, file=paste0(out,"/regions.bed"), row.names=F, col.names=F, quote=F, sep="\t" )
	regions<-paste0(out,"/regions.bed")

	#2- select probes with percentage of NA values > nalimit
	naprobes<-CpGNAexcl( betas,nalimit )
	filteredNAprobes<-length(naprobes)
	write.table(paste(naprobes, collapse="\n"), file=paste0(out, "/removed.txt"), row.names=F, sep="\t")

	#3- remove selected probes
	betas<-betas[ ! rownames(betas) %in% naprobes,]

	#4- remove remaining missing values
	print("Remove missing values")
	if ( sum(is.na(betas)) > 0 ) {

		if (missing=="mean"){
			betas<-replaceByMean(betas,pdata[,groups[1]])
		}
	
		if (missing=="impute"){
			betas<-imputeNA(betas,nalimit)
		}
	}

	#5- QC after processing
	print("QC after processing")
	for (group in groups) {

		jpeg(paste(out, group, "densityPlot2.jpg", sep="/"), width=800, height=800)
		densityPlot(betas, sampGroups = pdata[,group])
		dev.off()

		jpeg(paste(out, group, "densityBeanPlot2.jpg", sep="/"), width=800, height=800)
		densityBeanPlot(betas, sampGroups = pdata[,group])
		dev.off()

		jpeg(paste(out, group, "mdsPlot2.jpg", sep="/"), width=800, height=800)
		if (nrow(pdata) < 200 ){
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group],
						legendPos = "bottomleft",sampNames=samples)
		}else{
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group])
		}
		dev.off()
	}
	nbprobes2=nrow(betas)

	return( list(betas=betas, pdata=pdata, platform=platform, nbprobes1=nbprobes1, nbprobes2=nbprobes2, filteredFromList=filteredFromList, filteredNAprobes=c(), regions="") )

}
