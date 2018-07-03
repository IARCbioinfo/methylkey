#####################################################################################
#
# Title  : methylkey.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################

suppressWarnings(suppressPackageStartupMessages(library(GetoptLong))) 
suppressPackageStartupMessages(library(wateRmelon))
suppressPackageStartupMessages(library(data.table))
#library(RCurl)

#suppressWarnings(suppressPackageStartupMessages(library(methylumi)))
#suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))

#######################
#0- parsing options and inputs
nbc=6
idat="idat_rep"
package="minfi"
normalize="funnorm"
nalimit=0.2
missing="keep"
cell="None"
platform="unknown"
missing="mean"
GetoptLong(matrix(c(	"pdata=s", 	"pdata file",
			"idat=s",	"idat repository",
			"html=s",  	"report file, mandatory option",
			"out=s",  	"output dir",
			"samples=i",	"sample column index",
			"groups=i@",   	"group column index",
			"barcode=i",	"barcode column index",
			"cell=s",	"houseman or refFreeWasher",
			"nbc=i",	"number of cell type",
			"filter=s@",	"list of probes to filter",
			"normalize=s",	"Normalization method c(funnorm)",
			"nalimit=f",	"NA percentage cutoff to remove a probe",
			"missing=s",	"how to deal with missing values (NA) c(keep,mean,impute)",
			"package=s",	"minfi (default) or methylumi",
			"platform=s",	"illumina450k, illuminaEPIC, or matrix"
		), ncol=2, byrow=TRUE))

path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
print(path)
source( paste0(path, "/utils.r")  )
source( paste0(path, "/pca.r")  )
dir.create(out)

#######################
#1- Loading data
print("loading data")
pdata<-read.table(pdata, sep="\t", header=T)
pdata$sentrix_id<-as.factor(pdata$sentrix_id)
for (variable in colnames(pdata)){ if (length(levels(factor(pdata[,variable]))) < 4) { pdata[,variable]<-as.factor(pdata[,variable]) }  }
groups<-colnames(pdata)[groups]
samples<-pdata[,samples]

if (package=="minfi"){ 
	source(paste0(path,"/minfi.r"));
	pdata[,barcode]<-as.character(pdata[,barcode]) 
	colnames(pdata)[barcode]<-"Basename" 
}
if (package=="methylumi"){
 	source(paste0(path,"/methylumi.r"))
}

#readmeth pipeline change depending the loaded package
set<-readmeth(pdata=pdata, idat=idat, samples=samples, groups=groups, normalize=normalize, out=out, cell=cell)
print(set$platform)
print(platform)
if (set$platform != platform) { warning(paste0("Detected platform (" , set$platform, ") differ from expected (", platform, ") !")) }
platform=set$platform
pdata=set$pdata
betas=set$betas
nbprobes1=nrow(betas)

######################
#2- remove probes
print("Filter probes")
probes=c()
for (file in filter){
	probes <- unique( c(probes, fread(file)[[1]] ))
}
filteredFromList<-length(probes)

#select probes with percentage of NA values > CpGlimit
naprobes<-CpGNAexcl(betas,nalimit)
probes<-unique( c(probes, naprobes))
filteredNAprobes<-length(naprobes)

#remove select probes
totalFilteredProbes<-length(probes)
betas<-betas[ ! rownames(betas) %in% probes,]
nbprobes2<-nrow(betas)

######################
#3- remove missing values
print("Remove missing values")
if ( sum(is.na(betas)) > 0 ) {

	#Set remainings missing values to mean of betas for the probe.
	#@Novoloacaa
	if (missing=="mean"){
		probNAcont<-which(apply(betas,1,function(i) sum(is.na(i)))>0)
		for (i in names(probNAcont)){
			betas[i,is.na(betas[i,])]<-mean(betas[i,],na.rm=T)
		}
	}

	if (missing=="impute"){
		imputeNA(betas,nalimit)
	}
}
nbprobes3<-nrow(betas)

#######################
#4-PCA
print("running PCA")
suppressWarnings(  
tab_qvalue<-makepca( betas, pdata, out, colnames(pdata), 10 )
)

#######################
#5-DeltaBetas
print("Calculate Delta Betas")
save.image(file=paste0(out, "/meth.rdata") )
deltab<-getDeltaBetas(betas,pdata[,groups[1]])

#######################
#violin plot
violin_plot(pdata=pdata, betas=betas, group=groups[1], barcode=barcode, path=path, out=out, genome="hg19", cpg_annotated= NULL)

#######################
#5-Save
print("Save results")
write.table(pdata, file=paste0(out, "/pdata.txt"))
write.table(betas, file=paste0(out, "/betas.txt"), row.names=T, sep="\t")
write.table(paste(naprobes, collapse="\n"), file=paste0(out, "/removed.txt"), row.names=F, sep="\t")
write.table(deltab, file=paste0(out, "/deltabetas.txt"), row.names=T, sep="\t")
#save(betas, pdata, groups, samples, out, platform,  file=paste0(out, "/meth.rdata") )
save.image(file=paste0(out, "/meth.rdata") )

########################
#6-Create html report
print("Create report")

options<-paste(sapply(groups, function(x) sprintf("<option value=%s>%s</option>", x, x )), collapse="")[1]
group<-groups[1]

library(templates)
tpl<-tmpl(paste(readLines(paste0(path, "/report.tpl.html")), collapse="\n"))
tpl<-tmplUpdate(tpl, samples=nrow(pdata))

cat(tpl, file = html)






