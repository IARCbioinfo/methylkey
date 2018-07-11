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
GetoptLong(matrix(c(	"meth=s", 	"meth.rdata file path",
			"correction=s", "batch correction (sva, ssva or combat) , mandatory option",
			"batch=i",	"batch variable, optional",
			"variables=i@",	"vector of variables, mandatory option",
			"what=s",	"beta if new batch correction or mval if chained correction"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
print(path)
source( paste0(path, "/utils.r") )
source( paste0(path, "/pca.r") )

######################
#1-read data
print("load data")
load(meth)
process_id="D032"
comment="Correction from filtered betas"
if (what=="mval"){
	comment=paste0( "Correction chained from ", process_id )
}
process_id=idmaker(1)

######################
#2-check and reformat variables.
variables<-colnames(pdata)[variables]
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
}


if(correction == "combat"){

	print("running combat")	
	require(sva)
	mval<-ComBat(dat=mval, batch=pdata[,batch], mod=design, par.prior=TRUE, prior.plots=TRUE)
	#if (!exists("sva")){sva=NULL}
	correction=paste0(correction, " on batch ", colnames(pdata)[batch])
}


#######################
#7-PCA
print("running PCA")
formula1<-paste(formula1, collapse="") #pca id
suppressWarnings(  
tab<-makepca( mval, pdata, out, colnames(pdata), 10, id=process_id )
)
qvalue<-tab$qval
pvalue<-tab$pval

######################
#8-Create html report
print("save report")
library(templates)

save(betas, mval, pdata, groups, samples, design, platform, genome, html, out, process_id, file=paste0(out, "/meth.rdata") )

report<-paste(readLines(html), collapse="\n")
tpl<-tmpl(paste(readLines(paste0(path, "/batchcorrection.html.tpl")), collapse="\n"))
tpl<-tmplUpdate(tpl, qvalue=qvalue)

cat(report, tpl, file = html)



