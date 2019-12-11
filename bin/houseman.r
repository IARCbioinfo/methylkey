#######################################################
#     Houseman : estimate blood cell count            #
#######################################################
#

houseman<-function(RGset=RGset, pdata=pdata, referencePlatform="IlluminaHumanMethylation450k", out){
	
	if( referencePlatform=="IlluminaHumanMethylation450k" ){ require(FlowSorted.Blood.450k) }
	if( referencePlatform=="IlluminaHumanMethylationEPIC" ){ require(FlowSorted.Blood.epic) }

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

