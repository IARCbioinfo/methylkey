#####################################################################################
#
# Title		: dmps.r
# Author	: CahaisV@iarc.fr, RenardC@iarc.fr
# Date		: 07/10/2014
# Last Update	: 05/07/2018
#
#####################################################################################

library(GetoptLong)
library(limma)
library(qqman)

#######################
#parsing options and inputs
batch="no"
regions=""
model="None"
fdr="0.05"
genome="hg19"
method="ls"
nbpr2disp=50
output=NULL
celllines<-c("K562","Nhlf","Hsmm","Gm12878","Hmec","Hepg2","Huvec","Nhek","H1hesc")

GetoptLong(matrix(c(	"meth=s", 	"meth.rdata file path",
			"model=s",	"model",
			"variables=s@",	"vector of variables, mandatory option",
			"case=s",	"level to compare",
			"control=s",	"level to compare",
			"fdr=f",	"fdr cutoff for toptable",
			"method=s",	"method : ls or robust",
			"nbpr2disp=i",	"number of dmps to plot for heatmap",
			"celllines=s@",  "cellines",
			"output=s",	"output dir"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
templates<-paste0(path , "/../templates/")
datadir<-paste0(path , "/../data/")
print(path)
source( paste0(path, "/utils.r")  )
source( paste0(path, "/plots.r")  )
source( paste0(path, "/annot.r")  )

case=gsub("__cn__","",case) #chomp endline when the class parameter is the last one 
control=gsub("__cn__","",control) #chomp endline when the class parameter is the last one 


################################################
#1-read data
print("load data")
load(meth)
process_id=idmaker(1)
if ( !is.null(output)  ) { out<-output; dir.create(out) ; html=paste0(out,"/index.html")}
dir.create(paste0(out,"/",process_id))

#############################################
#2-Prepare variables
variables<-autoformat(variables)
variables<-colnames(pdata[variables])
pdata[ ,variables[1] ]<-as.factor(pdata[ ,variables[1] ])
for (v in variables) {
	if (is.factor(v)){
		pdata[ ,variables[v] ]<-make.names(pdata[ ,variables[v] ])
	}
}

#############################################
#3-Order the levels : control become intercept and case is compared to control.
case<-make.names(case)
control<-make.names(control)
pdata[,variables[1]] <- relevel(pdata[,variables[1]], case)
pdata[,variables[1]] <- relevel(pdata[,variables[1]], control)

#############################################
#4-create model
formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
if (model!="None"){ formula1 <- as.formula(paste("~ " ,model, collapse="") ) }
	
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
toptable<-topTable(eb, adjust="BH", number=Inf, p=fdr, sort.by="P")
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

#get regions
if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
if(regions == "") { stop( "\r Oups ! no region to annotate.") }

#complete annotation
cpg_annotated<-mk_annotate(regions,genome=genome)
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
if (nrow(toptable)<nbpr2disp){ nbpr2disp=nrow(toptable) }

if (nbpr2disp>1) {
	cpgs<-toptable$cpg[1:nbpr2disp]
	symbols<-toptable$annot.symbol[1:nbpr2disp]
	labrow<-paste0( cpgs, "(", symbols, ")")
	jpeg(paste(out, process_id, "heatmap.jpg",sep="/"), width=1600, height=1600)
	aheatmap(betas[cpgs,], annCol=pdata[,groups], labRow=labrow)
	dev.off()
}

########################
#11- barplots

p<-mk_barplot(annotated_regions=dm_annotated, betafc=as.numeric(dm$logFC), what="cpgi" )
ggsave(file=paste(out, process_id, "cpgi_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300, units=)

p<-mk_barplot(annotated_regions=dm_annotated, betafc=as.numeric(dm$logFC), what="genes" )
ggsave(file=paste(out, process_id, "genes_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)

########################
#12- histone mark
p<-mk_hmplot(dm, out, process_id, celllines)
ggsave(file=paste(out, process_id, "chromHMM_barplot.jpg", sep="/"), plot=p, width=20, height=20, dpi = 300)
#save.image(file=paste0(out,"/debug.rdata"))

########################
#MAPINFO
print("create info file")

info<-unique(data.frame(cpg_annotated)[,c(1,2,6)])
info<-merge(info, eb$p.value, by.x="name", by.y="row.names")
#info<-info[,c(1,2,3,8)]
colnames(info)<-c("TargetID", "CHR", "MAPINFO", "Pval")
#rownames(info)<-seq(1:nrow(info))

########################
#12-Create html report
print("Create report")
library(templates)

formula1<-paste(formula1, collapse=" ")
cmtx<-paste(cmtx, collapse=",")

dmps<-tmpl(paste(readLines(paste0(templates, "/dmps.html.tpl")), collapse="\n"))
dmps<-tmplUpdate(dmps, process_id=process_id)

analysis[process_id]=dmps

scripts<-tmpl(paste(readLines(paste0(templates, "/scripts.js")), collapse="\n"))
style<-tmpl(paste(readLines(paste0(templates, "/style.css")), collapse="\n"))
main<-tmpl(paste(readLines(paste0(templates, "/main.html.tpl")), collapse="\n"))
tpl<-paste(analysis, collapse="\n")
tpl<-tmplUpdate(main, scripts=scripts,style=style,tpl=tpl)

cat(tpl, "\n", file = html)

save(betas, mval, pdata, groups, samples, design, platform, method, analysis, genome, regions, toptable, html, out, process_id, info, file=paste0(out, "/meth.rdata") )




