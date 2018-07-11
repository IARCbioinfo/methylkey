#####################################################################################
#
# Title		: dmps.r
# Author	: CahaisV@iarc.fr, RenardC@iarc.fr
# Date		: 07/10/2014
# Last Update	: 05/07/2018
#
#####################################################################################

#library(GetoptLong)
#library(xtable)
#library(RColorBrewer)
#library(ggplot2 )
#suppressPackageStartupMessages(library(gplots))
#suppressPackageStartupMessages(library(methylumi))
#suppressPackageStartupMessages(library(wateRmelon))
#suppressPackageStartupMessages(library(qqman))
#suppressPackageStartupMessages(library(ggbio))
#suppressPackageStartupMessages(library(data.table))

library(GetoptLong)
library(limma)
library(qqman)
#library(annotatr)

#######################
#parsing options and inputs
batch="no"
regions=""
model="None"
qvalue="0.05"
genome="hg19"
GetoptLong(matrix(c(	"meth=s", 	"meth.rdata file path",
			"model=s",	"model",
			"variables=i@",	"vector of variables, mandatory option",
			"level1=s",	"level to compare",
			"level2=s",	"level to compare",
			"qvalue=f",	"qvalue cutoff for toptable",
			"method=s",	"method : ls or robust"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
print(path)
source( paste0(path, "/utils.r")  )

level1=gsub("__cn__","",level1) #chomp endline when the class parameter is the last one 
level2=gsub("__cn__","",level2) #chomp endline when the class parameter is the last one 

process_id=idmaker(1)

################################################
#1-read data
print("load data")
load(meth)
genome="hg19"
dir.create(paste0(out,"/",process_id))

#############################################
#2-Prepare variables
variables<-colnames(pdata)[variables]
pdata[ ,variables[1] ]<-as.factor(pdata[ ,variables[1] ])
for (v in variables) {
	if (is.factor(v)){
		pdata[ ,variables[v] ]<-make.names(pdata[ ,variables[v] ])
	}
}

#############################################
#3-Order the levels : level2 become intercept and level 1 is compared to level 2.
level1<-make.names(level1)
level2<-make.names(level2)
pdata[,variables[1]] <- relevel(pdata[,variables[1]], level1)
pdata[,variables[1]] <- relevel(pdata[,variables[1]], level2)

#############################################
#4-create model
formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
if (model!="None"){ formula1 <- as.formula(paste("~ " ,custum, collapse="") ) }
	
design<-model.matrix(formula1,data=pdata)

colnames(design)<-make.names(colnames(design))
cmtx <- makeContrasts(colnames(design)[2], levels= colnames(design))

#############################################
#5-remove cpgs when sum of its betas is 0
#it can produce a BUG when it is the case
iszero<-which(rowSums(mval) == 0 )
if(length(iszero) > 0){
	mval<-mval[-which(rowSums(mval) == 0 ),]
}

#############################################
#6-run model
print("running lmfit")
fit<-lmFit(mval,design, pdata, ndups=1, method=method)
rownames(cmtx)<-colnames(fit)
fitContrasts=contrasts.fit(fit,cmtx)
eb=eBayes(fitContrasts)
lmPvals = eb$p.value
chisq <- qchisq(1-eb$p.value,1)
lambda<-median(chisq)/qchisq(0.5,1)
#qqplot
jpeg(paste(out, process_id, "qqplot.jpg", sep="/"))
suppressWarnings( qq(lmPvals,title=paste("QQ plot: ","lambda=",lambda,sep="")) )
dev.off()

#############################################
#7-get top dmps
toptable<-topTable(eb, adjust="BH", number=Inf, p=qvalue, sort.by="P")
nbdmps=nrow(toptable)
print( paste0("you get ", nbdmps, " DMPS"))
if (nbdmps == 0) {
	#save(betas, mval, ssva, pdata, groups, samples, design, platform, html, out, file=meth )
	stop( "\r sorry !")
}

#############################################
#8-delta betas
print("Calculating delta betas")
deltabetas<-getDeltaBetas(betas[rownames(toptable),],pdata[,variables[1]])
toptable$deltaBetas<-deltabetas[,1]

#############################################
#9-annotation
print("Running annnotation")
annotations=NULL
if (genome=="hg19"){ 
	load(paste0(path,"/annotatr.hg19.rdata")) 
} else { 
	annotations<-buildAnnot(genome=genome) 
}

#get regions
if (platform=="IlluminaHumanMethylation450k"){regions=paste(path,"illumina450k.bed",sep="/")}
if (platform=="IlluminaHumanMethylationEPIC"){regions=paste(path,"illuminaEpic.bed",sep="/")}
if(regions == "") { stop( "\r Oups ! no region to annotate.") }

#complete annotation
cpg_annotated<-annotate(regions,genome=genome)
dm_annotated<-cpg_annotated[ cpg_annotated$name %in% rownames(toptable)]
dm<-merge(toptable, dm_annotated, by.x="row.names", by.y="name")
colnames(dm)[1]<-"cpg"
colnames(dm)[9]<-"chr"
dm$DM_Status[ as.numeric(dm$logFC) > 0 ] = "hyper"
dm$DM_Status[ as.numeric(dm$logFC) < 0 ] = "hypo"
write.table(dm, file=paste(out, process_id, "annotated.txt",sep="/"), row.names=F, sep="\t")

#toptable
dm_annotated<-makeGRangesFromDataFrame(dm, keep.extra.columns=TRUE)
dm<-dm[order(dm$annot.symbol),]
toptable<-dm[!duplicated(dm$cpg),c(9:11,1,13,2:8,23)]
toptable<-toptable[order(toptable$P.Value),]
write.table(toptable, file=paste(out, process_id, "toptable.txt",sep="/"), row.names=F, sep="\t")

#############################################
#10- Plots

#circusplot
print("circus plot")
gbm<-makeGRangesFromDataFrame(toptable)
gbm$meth.diff<-toptable$deltaBetas*100
p<-circusplot(gbm, genome=genome)
jpeg(paste(out, process_id, "circleplot.jpg", sep="/") , width=800, height=800)
print(p)
dev.off()

#heatmap
print("heatmap")
require(NMF)
nbpr2disp=50
if (nrow(toptable)>nbpr2disp){

	cpgs<-toptable$cpg[1:nbpr2disp]
	symbols<-toptable$annot.symbol[1:nbpr2disp]
	labrow<-paste0( toptable[cpgs], "(", symbols, ")")
	jpeg(paste(out, process_id, "heatmap.jpg",sep="/"), width=1600, height=1600)
	aheatmap(betas[cpgs,], annCol=pdata[,groups], labCol=samples)
	dev.off()

}

########################
#11-Create html report
print("Create report")
library(templates)

save(betas, mval, ssva, pdata, groups, samples, design, platform, genome, toptable, dm_annotated, html, out, file=paste0(out, "/meth.rdata") )

report<-paste(readLines(html), collapse="\n")
tpl<-tmpl(paste(readLines(paste0(path, "/dmps.html.tpl")), collapse="\n"))
formula1<-paste(formula1, collapse=" ")
cmtx<-paste(cmtx, collapse=",")
tpl<-tmplUpdate(tpl, process_id=process_id)
cat(report, tpl, file = html)

