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

readmeth<-function(pdata, idat , samples, groups, out, filters, nalimit=0.2, normalize="funnorm", cell="None"){
	
	print("minfi pipeline")
	require(minfi)
	require(wateRmelon)
	require(IlluminaHumanMethylationEPICmanifest)
	require(IlluminaHumanMethylation450kmanifest)

	#1-read idat files
	RGset <- read.metharray.exp(base = idat, targets=pdata, force=TRUE)
	betas<-getBeta(RGset)
	colnames(betas)<-samples
	platform<-RGset@annotation["array"]
	nbprobes1=nrow(betas)

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
		
		print("Funnorm")
		
		rawbetas=betas
		RGset<-preprocessFunnorm(RGset, sex=estSex$predictedSex)
		betas=getBeta(RGset)
		colnames(betas)<-samples

		#After normalization NA values are replace by values close to 0. This restore the NA status.
		probNAcont<-which(apply(rawbetas,1,function(i) sum(is.na(i)))>0)
		for (i in names(probNAcont)){
			betas[i,is.na(rawbetas[i,])]<-NA
		}
	}

	#6- select probes from list
	print("Filter probes")
	probes=c()
	for (file in filters){
		if(file.exists(file)){
			probes <- unique( c(probes, fread(file)[[1]] ))
		}
	}
	filteredFromList<-length(probes)

	#7- select probes with percentage of NA values > nalimit
	naprobes<-CpGNAexcl( betas,nalimit )
	probes<-unique( c(probes, naprobes) )
	filteredNAprobes<-length(naprobes)
	write.table(paste(naprobes, collapse="\n"), file=paste0(out, "/removed.txt"), row.names=F, sep="\t")

	#8- remove selected probes
	totalFilteredProbes<-length(probes)
	betas<-betas[ ! rownames(betas) %in% probes,]

	#9- remove remaining missing values
	print("Remove missing values")
	if ( sum(is.na(betas)) > 0 ) {

		if (missing=="mean"){
			betas<-replaceByMean(betas,pdata[,groups[1]])
		}
	
		if (missing=="impute"){
			betas<-imputeNA(betas,nalimit)
		}
	}

	#10- QC after processing
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
	
	return( list(betas=betas, pdata=pdata, platform=platform, nbprobes1=nbprobes1, nbprobes2=nbprobes2, filteredFromList=filteredFromList, filteredNAprobes=filteredNAprobes, regions="") )

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
	write.csv(wbc$counts, file=paste0(out, "/count.csv"), row.names=T)
	write.csv(wbc$compTable, file=paste0(out, "/compTable.csv"), row.names=T)
	
	#Add the cell proportions to the pData
	pdata=cbind(pdata,wbc$counts)
	
	#plot
	col=ncol(pdata)
	test<-pdata[,c(groups,(col-5):col)]
	suppressMessages(test.m <- melt(test))
	cols <- colorRampPalette(brewer.pal(5,"Dark2"))(length(levels(test[,1])))
	groups<-test.m[,1]
	
	jpeg(paste0(out, "/housemanplot.jpg"))
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
