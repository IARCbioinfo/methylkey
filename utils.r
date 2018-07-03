#####################################################################################
#
# Title  : utils.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 17/01/2018
#
#####################################################################################



#remove probes with percentage of missing values > CpGlimit
#@Novoloacaa
CpGNAexcl<-function(betas,CpGlimit=0.2){
	
	NArow<-apply(betas,1,function(i) sum(is.na(i)))
	exclprob<-rownames(betas)[NArow>CpGlimit*ncol(betas)]
	return(exclprob)
}


#impute missing values
imputeNA<-function(betas,nalimit){

	require(pamr)
	#impute missing values
	foo<-list( x=as.matrix(betas), y=colnames(betas) )
	suppressWarnings(betas<-pamr.knnimpute(foo, rowmax=nalimit , colmax=0.95)$x)
	return(betas)
}




#load and annotate RRBS regions
#loadRRBSRegion<-function(filename, genome="hg38"){

#	suppressPackageStartupMessages(require(annotatr))
#	dm_regions=read_regions(filename, genome=genome) 
#	#build annotation for genome
#	listOfAnnotations<-builtin_annotations()
#	annots<-listOfAnnotations[grep(genome, listOfAnnotations)]
#	annotations = build_annotations(genome = genome, annotations = annots)
#	#annotate cpgs 
#	suppressWarnings(dm_annotated = annotate_regions(dm_regions, annotations=annotations))
#	return(dm_annotated)
#}


#plot a circus manhattan plot
circusplot<-function(ranges, genome){

	if (genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38) ; species=BSgenome.Hsapiens.UCSC.hg38}
	if (genome=="hg19"){ require(BSgenome.Hsapiens.UCSC.hg19) ; species=BSgenome.Hsapiens.UCSC.hg19}
	print(genome)
	print(species)
	chr.len = seqlengths(species)
	chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
	myIdeo <- GRanges(seqnames = names(chr.len), ranges = IRanges(start = 1, chr.len))
	seqlevels(myIdeo) = names(chr.len)
	seqlevels(ranges) = names(chr.len)
	seqlengths(myIdeo)<-myIdeo@ranges@width
	seqlengths(ranges)<-myIdeo@ranges@width

	g.po <-ranges[ranges$meth.diff < 0 ]
	g.per<-ranges[ranges$meth.diff > 0 ]
	if (length(g.po) > 0 & length(g.per) > 0){
		values(g.po)$id = "hypo"
		values(g.per)$id = "hyper"

		p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
		p <- p + layout_circle(g.po, geom = "point", size = 1, aes(x = midpoint, y = meth.diff, color = id), radius = 19, trackWidth = 20) +		scale_colour_manual(values = c("magenta", "green"))
		p <- p + layout_circle(g.per, geom = "point", size = 1, aes(x = midpoint, y = meth.diff, color = id), radius = 41, trackWidth = 20)
		p<- p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), vjust = 0, radius = 55, trackWidth = 7)
		return(p)
	}
}



#check if sample cheet is complete for selected variables
checkpdata<-function(pdata, variables){

	for (v in variables) { 
		del<-which(is.na(pdata[,v]))
		if (length(del>0)){
			stop( "Error : NA values for ", paste(colnames(pdata)[v]) , ". Complete your pdata or choose an other variable", sep=" ")
		}

	}
}


#check if sample cheet is complete for selected variables
readpdata<-function(pdata){

	pdata<-read.table(pdata, sep="\t", header=T)
	for (v in 1:ncol(pdata)) { 
		del<-which(is.na(pdata[,v]))
		if (length(del>0))  { stop( "Error : NA values for ", paste(colnames(pdata)[v]) , ". Complete your pdata for", v , sep=" ") }
		if (!is.numeric(pdata[,v])) { pdata[,v]<-as.factor(pdata[,v]) } 
		return(pdata)
	}
}



#get genome annotation
getTxdb<-function(genome)
{
	if      (grepl("hg",genome)) { species="Hsapiens"  }
	else if (grepl("mm",genome)) { species="Mmusculus" }
	else { return ("no reference") }
	
	suppressMessages(pkg<-paste("TxDb", species, "UCSC", genome ,"knownGene", sep="."))
	print(pkg)
	library(pkg, character.only = TRUE)
	txdb=eval(parse(text = pkg))
	return(txdb)
}



#######################
# multiplot
multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#######################
# violin plot
violin_plot<-function(pdata,betas,group,barcode, path, out="./",genome="hg19",cpg_annotated=NULL){

	require(ggplot2)
	require(annotatr)
	require(data.table)

	if (!is.data.table(pdata)){ pdata<-data.table(pdata) }
	colnames(pdata)[barcode]<-"barcode"

	if (cpg_annotated==NULL){
		load(paste0(path,"/annotatr.hg19.rdata"))
		if (platform=="IlluminaHumanMethylation450k"){regions=paste0(path,"/illumina450k.bed")}
		if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(path,"/illuminaEpic.bed")}

		cpg_regions <- read_regions(regions, genome=genome, format="BED")
		cpg_annotated <- annotate_regions(cpg_regions, annotations=annotations, minoverlap = 0L)
	}

	annot<-data.table(cbind(cpg_annotated$name, cpg_annotated$annot$symbol, cpg_annotated$annot$type))
	annot<-annot[!grepl("chromatin", V3) ]
	annot<-unique(annot)
	annot[, c("ref","type","annot") := tstrsplit(V3,"_") ]

	colnames(betas)<-gsub("X","",colnames(betas))
	b<-data.table(melt(as.matrix(betas)))
	b<-merge(pdata[,.( barcode, grp=get(group))],b,by.x="barcode",by.y="Var2")

	cpg<-merge(b, annot[type=="cpg"], by.x="Var1", by.y="V1")
	genes<-merge(b, annot[type=="genes"], by.x="Var1", by.y="V1",allow.cartesian=TRUE)
	p1<-ggplot(cpg, aes_string(x="grp", y="value", fill="grp")) + geom_violin() + theme(axis.text.x=element_text(angle=50,hjust=1, size=8)) + guides(fill =FALSE) +
        ylab("Methylation %") + xlab("")
	p2<-ggplot(cpg, aes_string(x="annot", y="value", fill="grp")) + geom_violin() + theme(axis.text.x=element_text(angle=50,hjust=1, size=8)) + guides(fill =FALSE) +
        ylab("Methylation %") + xlab("")
	p3<-ggplot(genes, aes_string(x="annot", y="value", fill="grp")) + geom_violin(scale="area") + theme(axis.text.x=element_text(angle=50,hjust=1, size=8)) +
        ylab("Methylation %") + xlab("")

	jpeg(paste(out, "violin1.jpg", sep="/"), width=1000, height=800)
	p1
	dev.off()

	jpeg(paste(out, "violin2.jpg", sep="/"), width=1000, height=800)
	p2
	dev.off()

	jpeg(paste(out, "violin3.jpg", sep="/"), width=1000, height=800)
	p3
	dev.off()


	#jpeg(paste(out, "violin.jpg", sep="/"), width=1000, height=800)
	#multiplot(p1,p2,p3,cols=3)
	#dev.off()
}

#######################
# get Delta Betas

getDeltaBetas<-function(betas,group){

	means<-by(t(betas), group, colMeans)

	delta=list()
	for( level1 in levels(group) ) {
		for( level2 in levels(group) ) {
			if (level1 != level2){
				delta[[paste0(level1,"vs",level2)]]<-means[[level2]]-means[[level1]]
			}
		}
	}
	return(do.call(cbind, delta))
}













