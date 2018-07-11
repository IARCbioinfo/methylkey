#####################################################################################
#
# Title  : minfi.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################




#######################################################
#        MINFI PIPELINE                               #
#######################################################
#

readmeth<-function(pdata, idat , samples, groups, out, normalize="funnorm", cell="None"){
	
	print("minfi pipeline")
	require(minfi)
	require(wateRmelon)
	require(IlluminaHumanMethylationEPICmanifest)
	require(IlluminaHumanMethylation450kmanifest)

	#1-read idat files
	RGset <- read.metharray.exp(base = idat, targets=pdata, force=TRUE)
	betas<-getBeta(RGset)
	platform<-RGset@annotation["array"]

	#2-check gender
	GMsetEx <- mapToGenome(RGset) 
	estSex <- getSex(GMsetEx)
	pdata$predictedSex <- estSex$predictedSex
	pdata$predictedage <- agep(betas)
	
	#3- QC before processing
	nsamp<-nrow(pdata)

	jpeg(paste0(out, "/boxplot_colour1.jpg"), width=800, height=800)
	ylab<-"log2 intensity of both methylated and unmethylated probes"
	par(xaxt='n')
	boxplot(log2(RGset@assays$data$Red+1), col = "red", boxwex = 0.25, at= 1:nsamp - 0.175, ylab=ylab, labels=samples, cex=0.5)
	boxplot(log2(RGset@assays$data$Green+1), col = "green", boxwex = 0.25, at= 1:nsamp + 0.175, axis=F , add=T, cex=0.5)
	par(xaxt='s')
	axis(1, at=1:nsamp, labels=samples, tick=TRUE, las=2, cex.axis=0.8)
	dev.off()

	MSet <- preprocessRaw(RGset)
	qc <- getQC(MSet)
	jpeg(paste0(out, "/plotQC.jpg"), width=800, height=800)
	plotQC(qc)
	dev.off()

	jpeg(paste0(out, "/methylated1.jpg"), width=800, height=800)
	boxplot(log2(RGset@assays$data$Red+1), main="methylated (Red channel)", col = "red", las=2, cex.axi=0.8, labels=samples)
	dev.off()

	jpeg(paste0(out, "/unmethylated1.jpg"), width=800, height=800)
	boxplot(log(RGset@assays$data$Green+1), main="unmethylated (Green channel)", col = "green", las=2, cex.axi=0.8, labels=samples)
	dev.off()

	for (group in groups) {

		jpeg(paste(out, group, "densityPlot1.jpg", sep="/"), width=800, height=800)
		densityPlot(betas, sampGroups =pdata[,group])
		dev.off()

		jpeg(paste(out, group, "densityBeanPlot1.jpg", sep="/"), width=800, height=800)
		densityBeanPlot(betas, sampGroups = pdata[,group])
		dev.off()

		jpeg(paste(out, group, "mdsPlot1.jpg", sep="/"), width=800, height=800)
		if (nrow(pdata) < 200 ){
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group], legendPos = "bottomleft")
		}else{
			mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group])
		}
		dev.off()
	
	}

	#4-Houseman
	if(cell=="houseman"){ pdata<-houseman(RGset, data, out) }
	
	#5- Normalize with funnorm
	if(normalize=="funnorm"){
		rawbetas=betas
		RGset<-preprocessFunnorm(RGset, sex=estSex$predictedSex)
		betas=getBeta(RGset)

		#After normalization NA values are replace by values close to 0. This restore the NA status.
		probNAcont<-which(apply(rawbetas,1,function(i) sum(is.na(i)))>0)
		for (i in names(probNAcont)){
			betas[i,is.na(rawbetas[i,])]<-NA
		}
	}

	return( list(betas=betas, pdata=pdata, platform=platform ) )

}


#######################################################
#     Houseman : estimate blood cell count            #
#######################################################
#

houseman<-function(RGset=RGset, pdata=pdata, out){
	
	library(FlowSorted.Blood.450k)

	jpeg(paste0(out, "/ebcc.jpg"))
	wbc<-estimateCellCounts(RGset, returnAll=T, meanPlot=T, sex=as.factor(pdata$predictedSex))
	dev.off()

	#save count files
	write.csv(wbc$counts, file=paste0(output_dir, "/count.csv"), row.names=T)
	write.csv(wbc$compTable, file=paste0(out, "/compTable.csv"), row.names=T)
	
	#Add the cell proportions to the pData
	pdata=cbind(pdata,wbc$counts)
	
	#plot
	col=ncol(pdata)
	test<-pdata[,c(groups,(col-5):col)]
	suppressMessages(test.m <- melt(test))
	cols <- colorRampPalette(brewer.pal(5,"Dark2"))(length(levels(test[,1])))
	groups<-test.m[,1]
	
	jpeg(paste0(output_dir, "/housemanplot.jpg"))
	print(ggplot(test.m, aes(x=groups, y=value, fill = groups)) +
					geom_boxplot(outlier.colour = "black", outlier.size = 3) +
					scale_fill_manual(values = cols) +
					xlab("") +
					ylab("estimated wbc proportion") +
					theme( axis.text.x = element_blank() ) +
					facet_wrap(~ variable, scales = "free") +
					ggtitle(paste("Houseman cell type estimates by", colnames(test.m)[1] ,sep=" ")))
	
	dev.off()
	return(pdata)
}
