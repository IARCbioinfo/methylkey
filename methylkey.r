#####################################################################################
#
# Title  : methylkey.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################

suppressPackageStartupMessages(library(GetoptLong))
suppressPackageStartupMessages(library(wateRmelon))
suppressPackageStartupMessages(library(data.table))
#library(RCurl)

#suppressWarnings(suppressPackageStartupMessages(library(methylumi)))
#suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))

#######################
#0- parsing options and inputs
nbc=6
idat="idat_rep"
normalize="funnorm"
nalimit=0.2
missing="keep"
cell="None"
platform="unknown"
missing="mean"
genome="hg19"
barcode=0
regions=""
filter=c()
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
			"pipeline=s",	"rnbeads, minfi or methylumi",
			"platform=s",	"IlluminaHumanMethylation450k, IlluminaHumanMethylationEPIC, or matrix",
			"genome=s",	"reference genome"
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

#convert variables to factor
for (variable in colnames(pdata)){ 
	if (length(levels(factor(pdata[,variable]))) <= nrow(pdata)/8) 
	{ 
		pdata[,variable]<-as.factor(pdata[,variable]) 
	}  
}

#check groups
groups<-colnames(pdata)[groups]

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

#check samples
samples<-pdata[,samples]

#loading
if (pipeline=="minfi"){ 
	source(paste0(path,"/minfi.r"));
	pdata[,barcode]<-as.character(pdata[,barcode]) 
	colnames(pdata)[barcode]<-"Basename" 
}
if (pipeline=="methylumi"){
 	source(paste0(path,"/methylumi.r"))
}
if (pipeline=="rnbeads"){
 	source(paste0(path,"/tabular.r"))
}

#readmeth pipeline change depending the loaded pipeline
set<-readmeth(pdata=pdata, idat=idat, samples=samples, groups=groups, normalize=normalize, out=out, cell=cell)

if (set$platform != platform) 
{ 
	warning(paste0("Detected platform (" , set$platform, ") differ from expected (", platform, ") !")) 
}
if (pipeline=="rnbeads") 
{ 
	pdata$barcode=as.character(samples)
	regions=set$regions
}
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
naprobes<-CpGNAexcl( betas,nalimit )
probes<-unique( c(probes, naprobes) )
filteredNAprobes<-length(naprobes)

#remove selected probes
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

		meanBetas<-function(betas){
			rmeans<-matrix( rowMeans(betas,na.rm=T), ncol=ncol(betas), nrow=nrow(betas) )
			betas[is.na(betas) ] <- rmeans[is.na(betas)]
			return( betas )
		}

		for(level in pdata[groups[1]]){
			sel<-pdata[,groups[1]]==level
			betas[probNAcont,sel]<-meanBetas(betas[probNAcont,sel])
		}

	}

	if (missing=="impute"){
		imputeNA(betas,nalimit)
	}

}
nbprobes3<-nrow(betas)

#######################
#4- QC after processing
print("QC after processing")
for (group in groups) {

	jpeg(paste(out, group, "densityPlot2.jpg", sep="/"), width=800, height=800)
	densityPlot(betas, sampGroups = pdata[,group])
	dev.off()

	jpeg(paste(out, group, "densityBeanPlot2.jpg", sep="/"), width=800, height=800)
	densityBeanPlot(betas, sampGroups = pdata[,group])
	dev.off()

	jpeg(paste(out, group, "mdsPlot2.jpg", sep="/"), width=800, height=800)
	if (nrow(pdata) < 200 ){
		mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group],
					legendPos = "bottomleft",sampNames=samples)
	}else{
		mdsPlot(betas,numPositions=1000,sampGroups=pdata[,group])
	}
	dev.off()
}

#######################
#5-PCA
print("running PCA")
suppressWarnings(  
tab<-makepca( betas, pdata, out, colnames(pdata), 10 )
)
qvalue<-tab$qval
pvalue<-tab$pval

#######################
#6-DeltaBetas
print("Calculate Delta Betas")
deltab<-getDeltaBetas(betas,pdata[,groups[1]])

#######################
#7-violin plot
print("Drawing violin plot")
violin_plot(pdata=pdata, betas=betas, group=groups[1], barcode=barcode,platform=platform, path=path, out=out, genome=genome, regions=regions)

#######################
#8-Save
print("Save results")
process_id=idmaker(1)
write.table(pdata, file=paste0(out, "/pdata.txt"))
write.table(betas, file=paste0(out, "/betas.txt"), row.names=T, sep="\t")
write.table(paste(naprobes, collapse="\n"), file=paste0(out, "/removed.txt"), row.names=F, sep="\t")
write.table(deltab, file=paste0(out, "/deltabetas.txt"), row.names=T, sep="\t")
save(betas, pdata, groups, samples, platform, genome, html, out, process_id, file=paste0(out, "/meth.rdata") )
#save.image(file=paste0(out, "/meth.rdata") )

########################
#9-Create html report
print("Create report")
library(templates)
options<-paste(sapply(groups, function(x) sprintf("<option value=%s>%s</option>", x, x )), collapse="")[1]
group<-groups[1]

tpl<-tmpl(paste(readLines(paste0(path, "/report.html.tpl")), collapse="\n"))
tpl<-tmplUpdate(tpl, samples=nrow(pdata))

cat(tpl, "\n", file = html)






