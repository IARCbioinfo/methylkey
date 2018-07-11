#####################################################################################
#
# Title  : 450k-dmr.r
# Author : VargasH@iarc.fr, CahaisV@iarc.fr
# Date   : 22/03/2016
#
#####################################################################################

suppressPackageStartupMessages(library(GetoptLong))
#suppressWarnings(suppressPackageStartupMessages(library(wateRmelon)))
suppressPackageStartupMessages(library(DMRcate))
suppressPackageStartupMessages(library(DMRcatedata))
#suppressWarnings(suppressPackageStartupMessages(library(FDb.InfiniumMethylation.hg19)))
#suppressWarnings(suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)))
suppressPackageStartupMessages(library(RColorBrewer))
#suppressWarnings(suppressPackageStartupMessages(library(rtracklayer)))
#suppressWarnings(suppressPackageStartupMessages(library(ChIPseeker)))
suppressPackageStartupMessages(library(annotatr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
#data(dmrcatedata)


library(dmrcate.r)

#######################
#0- parsing options and inputs
GetoptLong(matrix(c(	"input=s",	"dmps.rdata",
			"analysis=s",	"type of analysis",
			"metho=s",	"method : ls or robust",
			"fdr=f",	"fdr cutoff",
			"pcutoff=f",	"pvalue cutoff",
			"report=s",	"html report",
			"outdir=s",	"report folder",
			"path=s",	"path to script files",
			"max=i",	"max"
		), ncol=2, byrow=TRUE))

source(paste(path,"utils.r",sep="/"))
dir.create(outdir)

########################
#1- read data
report_dir <- paste(sub(".dat$","",input), "files", sep="_")
dmps<-paste(report_dir, "dmps.rdata",sep="/")
load(dmps)
if (grep("Intercept" , colnames(design)[1] )) { colnames(design)[1] <- "(Intercept)" }

########################
#2- dmrcate 
if (platform!="sequencing"){

	genome="hg19"
	arannot=NULL
	niter = 25 ; if (metho=="robust"){ niter=200 }

	if (platform=="EPIC"){ arannot=c(array = "IlluminaHumanMethylationEPIC", annotation ="ilm10b2.hg19") }
	if (platform=="450K"){ arannot=c(array = "IlluminaHumanMethylation450k", annotation ="ilm10b2.hg19") }

	#Give betas to cpg.annotate (what="Betas"). They are converted to Mvalues in the function
	suppressWarnings(suppressMessages(
		myannotation<-cpg.annotate(object=betas,analysis.type=analysis, design=design, coef=2, datatype="array", what="Beta",fdr = fdr, annotation=arannot, arraytype=platform, method=metho, maxiter=niter)))
	suppressWarnings(suppressMessages(dmrcoutput<-dmrcate(myannotation,lambda=1000,C=2,pcutoff=pcutoff)))	
	suppressWarnings(suppressMessages(results.ranges <- extractRanges(dmrcoutput, genome = genome)))
}

if (platform=="sequencing"){
	
	final<-cbind(toptable$t, toptable$chr, toptable$start, toptable$logFC, toptable$adj.P.Val)
	colnames(final)<-c('stat', 'chr', 'pos', 'diff', 'fdr')
	suppressWarnings(suppressMessages(myannotation<-cpg.annotate( as.data.frame(final), datatype="sequencing", analysis.type="differential", method=metho, maxiter=niter)))	
	suppressWarnings(suppressMessages(dmrcoutput<-dmrcate(myannotation,lambda=1000,C=50,pcutoff=0.05)))
	dmrcoutput$results$coord<-paste("chr",dmrcoutput$results$coord,sep="")	
	dmrcoutput$results$coord<-sub("chr23", "chrX", dmrcoutput$results$coord )	
	dmrcoutput$results$coord<-sub("chr24", "chrY", dmrcoutput$results$coord )	
	suppressWarnings(suppressMessages(results.ranges <- extractRanges(dmrcoutput, genome = genome)))
}
nbdmrs<-length(results.ranges)
print( paste0("you get ", nbdmrs, " DMRS"))
if (nbdmrs == 0) stop( "\r sorry !")

#############################################
#3-annotation

n<-length(annotations[annotations$type==paste0(genome,"_cpg_shores")])
shores<-rep(c(paste0(genome,"_cpg_southShores"), paste0(genome,"_cpg_northShores")), (n+1)/2)
suppressWarnings(annotations[annotations$type==paste0(genome,"_cpg_shores")]$type<-shores[1:n])

n<-length(annotations[annotations$type==paste0(genome,"_cpg_shelves")])
shelves<-rep(c(paste0(genome,"_cpg_southShelves"), paste0(genome,"_cpg_northShelves")), (n+1)/2)
suppressWarnings(annotations[annotations$type==paste0(genome,"_cpg_shelves")]$type<-shelves[1:n])

suppressMessages(suppressWarnings(dmrs_annotated <- annotate_regions(results.ranges, annotations=annotations)))

#getCpGs
locs <- GRanges(myannotation$CHR, IRanges(myannotation$pos, myannotation$pos, names=myannotation$ID))
overlap <- findOverlaps(dmrs_annotated,locs)
overlap <- data.table( cbind( query=queryHits(overlap), subject=subjectHits(overlap) ) )
dmrs_annotated$cpg<-overlap[ , paste( names(locs)[subject], collapse="," ) , by=query ]$V1

foo<-as.numeric(dmrs_annotated$meanbetafc)
dmrs_annotated$DM_Status[ foo > 0 ] = "hyper"
dmrs_annotated$DM_Status[ foo < 0 ] = "hypo"
x_order = c('hyper','hypo')

#############################################
#4- outputs files
dmrs<-as.data.frame(dmrs_annotated, row.names=NULL)
colnames(dmrs)[1]<-"chr"
dmrs$dmrs<-paste(dmrs$chr, dmrs$start, dmrs$end, sep=":")

dup<-duplicated(dmrs$dmrs)
table<-dmrs[!dup,c(24,1:10,20,11,22)]

symbols<-dmrs[which(! is.na(dmrs$annot.symbol) ),c("dmrs","annot.symbol")]
symbols<-symbols[!duplicated(symbols$dmrs),]
rownames(symbols)<-symbols$dmrs
table$annot.symbol<-symbols[table$dmrs,"annot.symbol"]
	
write.table(dmrs, file=paste(outdir, "annotated.txt",sep="/"), row.names=F, sep="\t")
write.table(table, file=paste(outdir, "toptable.txt",sep="/"), row.names=F, sep="\t")
save.image(file=paste(outdir, "dmrs.rdata",sep="/"))

########################
#5- Create html report
print("Create html report")
write("<html><head></head><body>", file=report)
write("package = DMRCate <br/>", file=report, append=TRUE)
write(paste0("genome = ", genome, " <br/>"),file=report, append=TRUE)
write(paste0("analysis = ", analysis, " <br/>"),file=report, append=TRUE)
write(paste0("method = ", method, " <br/>"),file=report, append=TRUE)
write(paste0("fdr = ", fdr, " <br/>"),file=report, append=TRUE)
write(paste0("pcutoff = ", pcutoff, " <br/>"),file=report, append=TRUE)

write(paste("Differentially Methylated Regions : <a href='toptable.txt'>TopTable</a>"), file=report, append=TRUE)
write(paste(" (", nbdmrs  ,") "), file=report, append=TRUE)
write("<a href='annotated.txt'>Complete annotation</a><br/>", file=report, append=TRUE)
write("Continue to R : <a href='dmrs.rdata'>RDATA</a><br/>", file=report, append=TRUE)

write("<table><tr><td><a href='cpgIsland.jpg'><img src='cpgIsland.jpg' width='420' height='420'/></a></td>", file=report, append=TRUE)
write("<td><a href='genic.jpg'><img src='genic.jpg' width='420' height='420'/></a></td></tr>", file=report, append=TRUE)
write("<tr><td><a href='cpgcat1.jpg'><img src='cpgcat1.jpg' width='420' height='420'/></a></td>", file=report, append=TRUE)
write("<td><a href='genescat1.jpg'><img src='genescat1.jpg' width='420' height='420'/></a></td></tr>", file=report, append=TRUE)
write("<tr><td><a href='cpgcat2.jpg'><img src='cpgcat2.jpg' width='420' height='420'/></a></td>", file=report, append=TRUE)
write("<td><a href='genescat2.jpg'><img src='genescat2.jpg' width='420' height='420'/></a></td></tr></table><br/>", file=report, append=TRUE)

################################
#6- Annotation of DMRs
results.ranges@seqinfo<-annotations@seqinfo
#levels(seqnames(results.ranges))<-seqlevels(results.ranges)

suppressMessages( dm_annsum <- summarize_annotations(annotated_regions = dmrs_annotated) )
suppressMessages( dmrs_random_regions <- randomize_regions(regions = results.ranges) )
suppressMessages( dmrs_random_annotated <- annotate_regions(regions = dmrs_random_regions, annotations = annotations) )
suppressMessages( dmrs_annsum_rnd <- summarize_annotations(annotated_regions = dmrs_annotated, annotated_random = dmrs_random_annotated) )
print(dm_annsum)

################################
#6.1- cpgIsland annotation
annots_order = paste0( genome, c('_cpg_southShelves', '_cpg_southShores', '_cpg_islands', '_cpg_northShores', '_cpg_northShelves', '_cpg_inter'))
dm_vs_kg_annotations = plot_annotation(annotated_regions = dmrs_annotated, annotated_random = dmrs_random_annotated, x_label = 'CpG Island Annotations',y_label = 'Count', annotation_order = annots_order)
jpeg(paste(outdir, "cpgIsland.jpg",sep="/") , width=800, height=800)
print(dm_vs_kg_annotations)
dev.off()

#barplot v1 (default in annotatr)
dm_vs_cpg_cat1 = plot_categorical( annotated_regions = dmrs_annotated, x='DM_Status', fill='annot.type', x_order = x_order, fill_order = annots_order, position='fill', plot_title = 'DM Status by CpG Annotation Counts',legend_title = 'Annotations',x_label = 'DM status',y_label = 'Count', annotated_random = dmrs_random_annotated)
jpeg(paste(outdir, "cpgcat1.jpg",sep="/") , width=800, height=800)
print(dm_vs_cpg_cat1)
dev.off()

#barplot v2 (custom)
foo<-dmrs_annotated[dmrs_annotated$annot$type %in% annots_order,]
foo<-rbind(data.frame(type=as.factor(foo$annot$type), status=foo$DM_Status, val=1), data.frame(type=foo$annot$type, status="all", val=1) )
levels(foo$type)<-annots_order
colors<-brewer.pal(12,"Paired")[c(5,2,1,2,5,3)]
jpeg(paste(outdir, "cpgcat2.jpg",sep="/") , width=800, height=800)
ggplot(foo, aes(fill=type, x=status, y=val)) + geom_bar( stat="identity", position="fill") + coord_flip() + scale_fill_manual(values=colors)
dev.off()

################################
#6.2 Genic annotation
annots_order = paste0( genome, c('_genes_intergenic', '_genes_1to5kb', '_genes_5UTRs', '_genes_promoters', '_genes_firstexons', '_genes_exons', '_genes_introns'))
dm_vs_kg_annotations = plot_annotation(annotated_regions = dmrs_annotated, annotated_random = dmrs_random_annotated, x_label = 'Genics Annotations',y_label = 'Count', annotation_order = annots_order)
jpeg(paste(outdir, "genic.jpg",sep="/") , width=800, height=800)
print(dm_vs_kg_annotations)
dev.off()

#barplot v1 (default in annotatr)
dm_vs_cpg_cat2 = plot_categorical( annotated_regions = dmrs_annotated, x='DM_Status', fill='annot.type', x_order = x_order, fill_order = annots_order, position='fill', plot_title = 'DM Status by CpG Annotation Counts',legend_title = 'Annotations',x_label = 'DM status',y_label = 'Count', annotated_random = dmrs_random_annotated)
jpeg(paste(outdir, "genescat1.jpg",sep="/") , width=800, height=800)
print(dm_vs_cpg_cat2)
dev.off()

#barplot v2 (custom)
foo<-dmrs_annotated[dmrs_annotated$annot$type %in% annots_order,]
foo<-rbind(data.frame(type=as.factor(foo$annot$type), status=foo$DM_Status, val=1), data.frame(type=foo$annot$type, status="all", val=1) )
levels(foo$type)<-annots_order
colors<-brewer.pal(7, "Paired")
jpeg(paste(outdir, "genescat2.jpg",sep="/") , width=800, height=800)
ggplot(foo, aes(fill=type, x=status, y=val)) + geom_bar( stat="identity", position="fill") + coord_flip() + scale_fill_manual(values=colors)
dev.off()


########################
if (platform=="sequencing"){ write("<br/></body></html>", file=report, append=TRUE); stop( "\r") }
if (platform=="EPIC"){ suppressPackageStartupMessages(require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)) }
if (platform=="450K"){ suppressPackageStartupMessages(require(IlluminaHumanMethylation450kanno.ilmn12.hg19)) }
########################
#7- Plot dmrs
lvl<-levels(as.factor(pdata[,group]))
#color levels
x<-length(lvl)
myColors <- brewer.pal(x,"Set3")
names(myColors)<-lvl
myColors<-myColors[as.character(pdata[,group])]
#dmrs ranges
gdmr<-makeGRangesFromDataFrame(table)
colnames(betas)<-pdata[,samplenames]

#DMR.plot
for (dmr in 1:length(gdmr)){
	if (dmr > max) {break}
	name=paste(seqnames(gdmr[dmr]), gdmr[dmr]@ranges@start, gdmr[dmr]@ranges@start+gdmr[dmr]@ranges@width, sep="_")
	print(name)
	img<-paste(name, "dmrcate.jpg", sep=".")
	tryCatch({
		jpeg(paste(outdir, img ,sep="/"), width=800, height=800)
		suppressWarnings(DMR.plot(ranges=gdmr, dmr=dmr, CpGs=betas, phen.col=myColors, genome=genome, pch=20, toscale=T, plotmedians=T, arraytype=platform))
		dev.off()
		write(paste("<a href='", img, "'><img src='", img, "' width=400, height=400 /></a>", sep=""), file=report, append=TRUE)
		},error=function(cond) {print(paste(ref.cpg,cond,sep=" "))
	})
	
}

write("<br/></body></html>", file=report, append=TRUE)







