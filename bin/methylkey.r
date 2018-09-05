#!/usr/bin/Rscript
#####################################################################################
#
# Title  : methylkey_load.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################

suppressPackageStartupMessages(library(GetoptLong))

#######################
#0- parsing options and inputs
nbc=6
idat="idat_rep"
pipeline="minfi"
normalize="funnorm"
nalimit=0.2
cell="None"
platform="unknown"
missing="mean"
genome="hg19"
barcode="barcode"
out="methylkey"
regions=""
violin=FALSE
filters=c()

GetoptLong(matrix(c(	"pdata=s", 	"pdata file",
			"idat=s",	"idat repository",
			"out=s",  	"output file directory",
			"samples=s",	"sample column name or index",
			"groups=s@",   	"group column name or index",
			"platform=s",	"IlluminaHumanMethylation450k, IlluminaHumanMethylationEPIC, or matrix",
			"barcode=s",	"barcode column name or index for array data",
			"pipeline=s",	"minfi, methylumi or rnbeads",
			"cell=s",	"houseman or refFreeWasher",
			"nbc=i",	"number of cell type",
			"normalize=s",	"Normalization method c(funnorm)",
			"filters=s@",	"list of probes to filter",
			"nalimit=f",	"NA percentage cutoff to remove a probe",
			"missing=s",	"how to deal with missing values (keep,mean,impute)",
			"genome=s",	"reference genome",
			"violin",	"draw violin plots"
		), ncol=2, byrow=TRUE))


suppressPackageStartupMessages(library(wateRmelon))
suppressPackageStartupMessages(library(data.table))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2]) 
templates<-paste0(path , "/../templates/")
datadir<-paste0(path , "/../data/")

source( paste0(path, "/utils.r")  )
source( paste0(path, "/plots.r")  )
source( paste0(path, "/annot.r")  )
source( paste0(path, "/pca.r")  )
source( paste0(path, "/missingValues.r")  )

print(out)
dir.create(out)
html=paste0(out,"/index.html")

samples<-suppressWarnings( autoformat(samples) )
barcode<-suppressWarnings( autoformat(barcode) )
groups <-suppressWarnings( autoformat(groups) )
if( !is.null(filters) ){ filters<-unlist(strsplit(filters,",")) }

process_id=idmaker(1)
print(process_id)

#######################
#1- Loading data
print("loading data")
pdata<-read.table(pdata, sep="\t", header=T)

#convert variables to factor
for (variable in colnames(pdata)){ 
	if (length(levels(factor(pdata[,variable]))) <= nrow(pdata)/8) 
	{ 
		pdata[,variable]<-as.factor(pdata[,variable]) 
	}  
}

#check samples
samples<-pdata[samples]
samplesnames<-pdata[,colnames(samples)]

#check barcode
if(barcode != "barcode"){
	barcode<-colnames(pdata[barcode])
}

#check groups
groups<-colnames(pdata[groups])
print(groups)

for (group in groups) {
	dir.create(paste(out, group, sep="/"))
}

for (i in 1:length(groups)) 
{
	if ( !is.factor(pdata[,groups[i]] ))
	{
		groups<-groups[-i]
		print(paste0("warning : ", groups[i], " removed from groups because it is numeric." ))
	}
}

#######################
#2- Running pipeline
if (pipeline=="minfi"){ 
	source(paste0(path,"/pipeline_minfi.r"));
	pdata[,barcode]<-as.character(pdata[,barcode]) 
	pdata$Basename<-pdata[,barcode] 
	tmplfile<-paste0(templates, "/pipeline_minfi.html.tpl")
}
if (pipeline=="methylumi"){
 	source(paste0(path,"/pipeline_methylumi.r"))
	tmplfile<-paste0(templates, "/pipeline_methylumi.html.tpl")
	if(normalize=="funnorm"){normalize="bmiq"} # default for methylumi is bmiq
}
if (pipeline=="rnbeads"){
 	source(paste0(path,"/pipeline_rnbeads.r"))
	tmplfile<-paste0(templates, "/pipeline_rnbeads.html.tpl")
	pdata$barcode=as.character(samples)
}

#readmeth pipeline change depending the loaded pipeline
data<-suppressPackageStartupMessages( readmeth(pdata=pdata, idat=idat, samples=samplesnames, groups=groups, normalize=normalize, filters=filters, nalimit=nalimit, out=out, cell=cell) )

regions=data$regions
platform=data$platform
pdata=data$pdata
betas=data$betas
nbprobes1=data$nbprobes1
nbprobes2=data$nbprobes2
filteredNAprobes=data$filteredNAprobes
filteredFromList=data$filteredFromList

#######################
#3- PCA
print("running PCA")
pca<-suppressPackageStartupMessages( makepca( betas, pdata, out, colnames(pdata), nPC=10, id=process_id ) )
pvalue<-pca$pvalue
qvalue<-pca$qvalue

#######################
#4- DeltaBetas
#save.image(file=paste0(out, "/debug.rdata") )
print("Calculate Delta Betas")
deltab<-""
deltab<-getDeltaBetas(betas,pdata[,groups[1]])

#######################
#5-violin plot
if (violin){
	print("Drawing violin plot")
	save.image(file=paste0(out, "/debug.rdata") )
	violin_plot(pdata=pdata, betas=betas, group=groups[1], samples=colnames(samples), platform=platform, path=path, out=out, genome=genome, regions=regions)
}
#######################
#6- Save
print("Save results")
write.table(pdata, file=paste0(out, "/pdata.txt"))
write.table(betas, file=paste0(out, "/betas.txt"), row.names=T, sep="\t")
write.table(deltab, file=paste0(out, "/deltabetas.txt"), row.names=T, sep="\t")
#save.image(file=paste0(out, "/debug.rdata") )

########################
#7- Create html report
print("Create report")
library(templates)

if (is.null(filters)){ filter="none"} else {filters<-paste0(basename(filters),collapse=", ")}
correction = "on filtered/normalized betas"
options<-paste(sapply(groups, function(x) sprintf("<option value=%s>%s</option>", x, x )), collapse="")[1]
tpl<-tmpl(paste(readLines(tmplfile), collapse="\n"))
tpl<-tmplUpdate( tpl, samples=nrow(pdata), group=groups[1] )

if (violin){
	vio<-tmpl(paste(readLines(paste0(templates, "/violin.html.tpl")), collapse="\n"))
	vio<-tmplUpdate( vio, out=out )
	tpl<-paste(tpl,vio,collapse="\n") 
}

pca<-tmpl(paste(readLines(paste0(templates, "/pca.html.tpl")), collapse="\n"))
pca<-tmplUpdate( pca, pvalue=pvalue, qvalue=qvalue )
tpl<-paste(tpl,pca,collapse="\n") 

analysis<-list()
analysis[process_id]=tpl

scripts<-tmpl(paste(readLines(paste0(templates, "/scripts.js")), collapse="\n"))
style<-tmpl(paste(readLines(paste0(templates, "/style.css")), collapse="\n"))
main<-tmpl(paste(readLines(paste0(templates, "/main.html.tpl")), collapse="\n"))
tpl<-tmplUpdate(main, scripts=scripts,style=style,tpl=tpl)

cat(tpl, "\n", file = html)

save(betas, pdata, groups, samples, platform, genome, regions, html, out, analysis, file=paste0(out, "/meth.rdata") )





