#!/usr/bin/Rscript
#####################################################################################
#
# Title  : batchcorrection.r
# Author : CahaisV@iarc.fr
# Date   : 04/07/2018
#
#####################################################################################

suppressPackageStartupMessages(library(GetoptLong))
suppressPackageStartupMessages(library(wateRmelon))

batch="no"
what="betas"
correction="sva"
out=NULL
GetoptLong(matrix(c(	"meth=s", 	"meth.rdata file path",
			"correction=s", "batch correction (sva, ssva or combat) , mandatory option",
			"batch=s",	"batch variable, optional",
			"variables=s@",	"vector of variables, mandatory option",
			"what=s",	"beta if new batch correction or mval if chained correction",
			"out=s",	"out dir"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
templates<-paste0(path , "/../templates/")
print(path)
source( paste0(path, "/utils.r") )
source( paste0(path, "/pca.r") )

######################
#1-read data
print("load data")
load(meth)
comment="Correction from filtered betas"
if (what=="mval"){
	comment=paste0( "Correction chained from ", process_id )
}
process_id=idmaker(1)
if ( !is.null(out)  ) { out<-out; dir.create(out) ; html=paste0(out,"/index.html")}
suppressWarnings( batch<-autoformat(batch) )

######################
#2-check and reformat variables.
variables<-autoformat(variables)
variables<-colnames(pdata[variables])
for (v in variables) {
	pdata[ ,v ]<-gsub(" ",".",pdata[ ,v ] )
	pdata[ ,v ]<-gsub("-","_",pdata[ ,v ] )
}

######################
#3-create model
formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
design<-model.matrix(formula1,data=pdata)

######################
#4-mvalues
if(what=="betas"){
	mval<-beta2m(as.matrix(betas))
	mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
}

######################
#5-Remove a probe when sum of betas is 0.
iszero<-which(rowSums(mval) == 0 )
if(length(iszero) > 0){
	mval<-mval[-which(rowSums(mval) == 0 ),]
}

######################
#6-Batch correction
if(correction == "sva"){

	print("running sva")

	require(sva)
	
	if (length(variables)==1){ 
		sva<-sva(mval,design)
	}
	else{
		formula0 <- as.formula(paste("~ 1+", paste(variables[-1],collapse="+")))
		model0<-model.matrix(formula0, data=pdata)
		sva<-sva(mval,design,model0)
	}
	#design<-cbind(design,sva$sv)
	mval = t(residuals(lm(t(mval)~sva$sv)))

	formula1<-paste(formula1, collapse="")
	correction=paste0("on SVA corrected mvalues with model ", formula1)
}

if(correction == "ssva"){
	
	print("running SmartSVA")

	require(SmartSVA)

	#estimate number of surrogates
	formula <- as.formula(paste("t(betas) ~ " ,paste(variables,collapse="+")))
	Y.r <- t(resid(lm(formula, data=pdata)))
	n.sv <- EstDimRMT(Y.r, FALSE)$dim
	
	if (length(variables)==1){ 
		ssva<-smartsva.cpp(mval, design, n.sv=n.sv)
	}
	else{
		formula0 <- as.formula(paste("~ 1+", paste(variables[-1],collapse="+")))
		model0<-model.matrix(formula0, data=pdata)
		ssva<-smartsva.cpp(mval, design, model0, n.sv=n.sv)
	}
	#design<-cbind(design,sva$sv)
	mval = t(residuals(lm(t(mval)~sva$sv)))
	
	formula1<-paste(formula1, collapse="")
	correction=paste0("on SSVA corrected mvalues with model ", formula1)
}


if(correction == "combat"){

	print("running combat")	
	require(sva)
	mval<-ComBat(dat=mval, batch=pdata[,batch], mod=design, par.prior=TRUE, prior.plots=TRUE)
	#if (!exists("sva")){sva=NULL}
	correction=paste0("on combat corrected mvalues for batch ", colnames(pdata)[batch] )
}


#######################
#7-PCA
print("running PCA")
 
tab<-makepca( mval, pdata, out, colnames(pdata), nPC=10, id=process_id ) 
qvalue<-tab$qval
pvalue<-tab$pval

######################
#8-Create html report
print("save report")
library(templates)

batch<-tmpl(paste(readLines(paste0(templates, "/pca.html.tpl")), collapse="\n"))
batch<-tmplUpdate(batch, qvalue=qvalue)
analysis[process_id]=batch

scripts<-tmpl(paste(readLines(paste0(templates, "/scripts.js")), collapse="\n"))
style<-tmpl(paste(readLines(paste0(templates, "/style.css")), collapse="\n"))
main<-tmpl(paste(readLines(paste0(templates, "/main.html.tpl")), collapse="\n"))
tpl<-paste(analysis, collapse="\n")
tpl<-tmplUpdate(main, scripts=scripts,style=style,tpl=tpl)

cat(tpl, "\n", file = html)

save(betas, mval, pdata, groups, samples, design, platform, genome, regions, analysis, html, out, process_id, file=paste0(out, "/meth.rdata") )






