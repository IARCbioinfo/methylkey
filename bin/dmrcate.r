#!/usr/bin/Rscript
#####################################################################################
#
# Title  : dmrcate.r
# Author : CahaisV@iarc.fr
# Last update : 19/07/2018
#
#####################################################################################

library(GetoptLong)
library(limma)
library(DMRcate)
library(data.table)

#######################
#0- parsing options and inputs
type="differential"
fdr=0.05
pcutoff=1e-5
niter=25
output=NULL
GetoptLong(matrix(c(	"meth=s",	"meth.rdata file path",
			"type=s",	"type of analysis",
			"fdr=f",	"fdr cutoff",
			"pcutoff=f",	"pvalue cutoff",
			"niter=i",	"number of iteration",
			"output=s",	"output dir"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
templates<-paste0(path , "/../templates/")
datadir<-paste0(path , "/../data/")
print(path)
source( paste0(path, "/utils.r")  )
source( paste0(path, "/plots.r")  )
source( paste0(path, "/annot.r")  )

################################################
#1-read data
print("load data")
load(meth)
if ( !is.null(output)  ) { out<-output; dir.create(out) ; html=paste0(out,"/index.html")}
process_id=idmaker(1)
dir.create(paste0(out,"/",process_id))

cmtx <- makeContrasts(colnames(design)[2], levels= colnames(design))
colnames(design)[1] <- "(Intercept)"

if (method=="robust"){ niter=200 }

########################
#2- dmrcate 
if (platform=="IlluminaHumanMethylationEPIC" | platform=="IlluminaHumanMethylation450k"){

	arannot=NULL
	if (platform=="IlluminaHumanMethylationEPIC"){ arraytype="EPIC"; arannot=c(array = "IlluminaHumanMethylationEPIC", annotation ="ilm10b2.hg19") }
	if (platform=="IlluminaHumanMethylation450k"){ arraytype="450K"; arannot=c(array = "IlluminaHumanMethylation450k", annotation ="ilm10b2.hg19") }

	myannotation<-cpg.annotate(object=mval, analysis.type=type, design=design, coef=2, datatype="array", arraytype=arraytype, what="M", fdr=fdr, annotation=arannot,  method=method, maxiter=niter, cont.matrix=cmtx)	

	dmrcoutput<-dmrcate(myannotation,lambda=1000,C=2,pcutoff=pcutoff)	
	results.ranges <- extractRanges(dmrcoutput, genome=genome)
}

if (platform=="sequencing"){
	
	final<-cbind(toptable$t, toptable$chr, toptable$start, toptable$logFC, toptable$adj.P.Val)
	colnames(final)<-c('stat', 'chr', 'pos', 'diff', 'fdr')
	suppressWarnings(suppressMessages(myannotation<-cpg.annotate( as.data.frame(final), datatype="sequencing", analysis.type="differential", method=method, maxiter=niter)))	
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
print("Running annnotation")

dmrs_annotated<-mk_annotate(results.ranges,genome=genome)

#overlapping.probes
locs <- GRanges(myannotation$CHR, IRanges(myannotation$pos, myannotation$pos, names=myannotation$ID))
overlap <- findOverlaps(dmrs_annotated,locs)
overlap <- data.table( cbind( query=queryHits(overlap), subject=subjectHits(overlap) ) )
dmrs_annotated$overlapping.probes<-overlap[ , paste( names(locs)[subject], collapse="," ) , by=query ]$V1

#dmrs
dmrs<-data.frame(dmrs_annotated, row.names=NULL)
colnames(dmrs)[1]<-"chr"
dmrs$dmrs<-paste(dmrs$chr, dmrs$start, dmrs$end, sep=":")

#remove duplicated lines
dmrs<-dmrs[order(dmrs$annot.symbol),]
dup<-duplicated(dmrs$dmrs)
table<-dmrs[!dup,c(23,1:10,20,11,22)]

#save files	
dmrs<-dmrs[order(dmrs$minfdr),]
table<-table[order(table$minfdr),]
write.table(dmrs, file=paste(out, process_id, "annotated.txt",sep="/"), row.names=F, sep="\t")
write.table(table, file=paste(out, process_id, "toptable.txt",sep="/"), row.names=F, sep="\t")


save.image(file=paste0(out, "/debug.rdata") )


########################
#barplots
save.image(file=paste0(out, "/debug.rdata"))

p<-mk_barplot(annotated_regions=data.frame(dmrs_annotated), betafc=as.numeric(dmrs_annotated$meanbetafc), what="cpgi" )
ggsave(file=paste(out, process_id, "cpgi_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300, units=)

p<-mk_barplot(annotated_regions=data.frame(dmrs_annotated), betafc=as.numeric(dmrs_annotated$meanbetafc), what="genes" )
ggsave(file=paste(out, process_id, "genes_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)


########################
#11-Create html report
print("Create report")
library(templates)

formula1<-paste(colnames(design)[-1], collapse=" + ")
cmtx<-paste(cmtx, collapse=",")

dmrs<-tmpl(paste(readLines(paste0(templates, "/dmrcate.html.tpl")), collapse="\n"))
dmrs<-tmplUpdate(dmrs, process_id=process_id)

analysis[process_id]=dmrs

scripts<-tmpl(paste(readLines(paste0(templates, "/scripts.js")), collapse="\n"))
style<-tmpl(paste(readLines(paste0(templates, "/style.css")), collapse="\n"))
main<-tmpl(paste(readLines(paste0(templates, "/main.html.tpl")), collapse="\n"))
tpl<-paste(analysis, collapse="\n")
tpl<-tmplUpdate(main, scripts=scripts,style=style,tpl=tpl)

cat(tpl, "\n", file = html)

save(betas, mval, pdata, groups, samples, design, platform, method, analysis, genome, regions, toptable, html, out, process_id, info, table, file=paste0(out, "/meth.rdata") )













