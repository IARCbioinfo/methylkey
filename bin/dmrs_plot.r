#!/usr/bin/Rscript
#####################################################################################
#
# Title		: dmps.r
# Author	: CahaisV@iarc.fr
# Date		: 24/07/2018
# Last Update	: 24/07/2018
#
#####################################################################################

library(GetoptLong)
library(RColorBrewer)
suppressPackageStartupMessages(library(GenomicRanges))


dmrs=c()
max=0
output=NULL
#######################
#parsing options and inputs
GetoptLong(matrix(c(	"meth=s", 	"meth.rdata file path",
			"max=i",	"max dmr to plot",
			"dmrs=s@",	"dmrs to plot",
			"output=s",	"output dir"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
templates<-paste0(path , "/../templates/")
print(path)

if ( max==0 & length(dmrs)==0) { stop("no dmr to plot") }

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
group<-groups[1]

########################
# Plot DMRs function

mkPlotDMR<-function(dmr){

	require(DMRcate)
	dmrplot=""
	dmrname=table$dmrs[dmr]
	print(dmrname)
	img<-paste0(process_id, "/", dmrname, ".dmrplot.jpg")
	tryCatch({
		jpeg(paste(out, img ,sep="/"), width=800, height=800)
		DMR.plot(ranges=gdmr, dmr=dmr, CpGs=betas, phen.col=myColors, genome=genome, pch=20, toscale=T, plotmedians=T, arraytype=arraytype)
		dev.off()
		dmrplot<-paste("<a href='", img, "'><img src='", img, "' width=400, height=400 /></a>", sep="")
	},error=function(cond) { dmrplot<-paste(dmrname, cond, sep=" ") })
	return(dmrplot)

}

########################
#

if (platform=="sequencing"){ }
if (platform=="IlluminaHumanMethylationEPIC"){ arraytype="EPIC"; require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) }
if (platform=="IlluminaHumanMethylation450k"){ arraytype="450K"; require(IlluminaHumanMethylation450kanno.ilmn12.hg19) }

#######################
#3- colors
lvl<-levels(as.factor(pdata[,group]))
myColors <- brewer.pal(length(lvl),"Set1")
names(myColors)<-lvl
myColors<-myColors[as.character(pdata[,group])]

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

if (max > length(toptable$cpg) ) { max=length(toptable$cpg) }
for (dmr in 1:max){
	
	dmrplots[dmr]<-mkPlotDMR(dmr)	
}

dmrplots<-paste(dmrplots, collapse="\n")

########################
#3-Create html report
print("Create report")
library(templates)

dmrs<-tmpl(paste(readLines(paste0(templates, "/dmrplot.html.tpl")), collapse="\n"))
dmrs<-tmplUpdate(dmrs, process_id=process_id)

analysis[process_id]=dmrs

scripts<-tmpl(paste(readLines(paste0(templates, "/scripts.js")), collapse="\n"))
style<-tmpl(paste(readLines(paste0(templates, "/style.css")), collapse="\n"))
main<-tmpl(paste(readLines(paste0(templates, "/main.html.tpl")), collapse="\n"))
tpl<-paste(analysis, collapse="\n")
tpl<-tmplUpdate(main, scripts=scripts,style=style,tpl=tpl)

cat(tpl, "\n", file = html)

save(betas, mval, pdata, groups, samples, design, platform, method, analysis, genome, regions, toptable, html, out, process_id, info, table, file=paste0(out, "/meth.rdata") )







