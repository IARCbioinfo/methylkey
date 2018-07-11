#####################################################################################
#
# Title  : methylumi.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 04/07/2018
#
#####################################################################################


#######################################################
#        METHYLUMI PIPELINE                           #
#######################################################
#

readmeth<-function(pdata, idat , samples, groups, out, normalize="bmiq", cell="None"){

	print("methylumi pipeline")
	require(methylumi)
	require(wateRmelon)
	require(RPMM)
	require(IlluminaHumanMethylationEPICmanifest)
	require(IlluminaHumanMethylation450kmanifest)
	mc.cores=4

	#1-read idat files
	methylumi<-readEPIC(pdat=pdata, idatPath=idat, parallel=T)
	betas<-betas(methylumi)
	sampleNames(methylumi)<-samples
	platform=gsub("Epic","EPIC",methylumi@annotation)
	
	#2-QC before processing
	jpeg(paste0(out, "/boxplot_colour1.jpg"), width=800, height=800)
	boxplotColorBias(methylumi)
	dev.off()

	jpeg(paste0(out, "/methylated1.jpg"), width=800, height=800)
	boxplot(log(methylated(methylumi)+1), col = "red", main="methylated (Red channel)", las=2, cex.axi=0.8)
	dev.off()

	jpeg(paste0(out, "/unmethylated1.jpg"), width=800, height=800)
	boxplot(log(unmethylated(methylumi)+1), col = "green", main="unmethylated (Green channel)", las=2, cex.axi=0.8)
	dev.off()

	for (group in groups) {

		jpeg(paste(out, group, "densityPlot1.jpg", sep="/"), width=800, height=800)
		densityPlot(betas, sampGroups = pdata[,group])
		dev.off()

		jpeg(paste(out, group, "/densityBeanPlot1.jpg", sep="/"), width=800, height=800)
		densityBeanPlot(betas, sampGroups = pdata[,group])
		dev.off()

		jpeg(paste(out, group, "/mdsPlot1.jpg", sep="/"), width=800, height=800)
		if (dim(pdata)[1] < 200 ){
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group], legendPos = "bottomleft",sampNames=pdata[,group])
		}else{
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group])
		}
		dev.off()
	}

	#3-normalize
	methylumi<-pfilter(methylumi, perc=100)
	if(normalize=="bmiq" ){
		methylumi<-BMIQ(methylumi)
	}
	methylumi<-adjColorBias.quantile(methylumi)
	betas<-betas(methylumi)
	
	return( list(betas=betas, pdata=pdata, platform=platform ) )
}

