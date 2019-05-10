#!/usr/bin/Rscript
#####################################################################################
#
# Title  : methylkey.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################

suppressPackageStartupMessages(library(GetoptLong))

#######################
#0- parsing options and inputs
pdata="pdata"
samples="Sample_Name"
groups="Sample_Group"
barcode="barcode"
idat="idat"
pipeline="minfi"
normalize="funnorm"
nalimit=0.2
missing="mean"
cell="None"
nbc=6
genome="hg19"
violin=FALSE
filters=c("default")
variables=c()
case=""
control=""
batch="no"
correction="sva"
meth="meth.rdata"
out="methylkey"
formula="None"
fdr="0.05"
genome="hg19"
method="ls"
hsize=50
celllines="default"
type="differential"
pcutoff=1e-5
niter=25

dmps=c()
dmrs=c()
max=0
win=10000
separator="\t"

GetoptLong(matrix(c("tool|T=s", "tool",	
      "pdata=s", 	"pdata file",
			"idat=s",	"idat repository",
			"out=s",  	"output file directory",
			"samples=s",	"sample column name or index",
			"groups=s@",   	"group column name or index",
			"barcode=s",	"barcode column name or index for array data",
			"pipeline=s",	"minfi, methylumi or rnbeads",
			"cell=s",	"houseman or refFreeWasher",
			"nbc=i",	"number of cell type",
			"normalize=s",	"Normalization method c(funnorm)",
			"filters=s@",	"list of probes to filter",
			"nalimit=f",	"NA percentage cutoff to remove a probe",
			"missing=s",	"how to deal with missing values (keep,mean,impute)",
			"genome=s",	"reference genome",
			"violin",	"draw violin plots",
			"meth=s", 	"meth.rdata file path",
			"correction=s", "batch correction (sva, ssva or combat) , mandatory option",
			"batch=s",	"batch variable, optional",
			"variables=s@",	"vector of variables, mandatory option",
			"formula=s",	"custom or advance model",
			"case=s",	"level to compare",
			"control=s",	"level to compare",
			"fdr=f",	"fdr cutoff for toptable",
			"method=s",	"method : ls or robust",
			"hsize=i",	"number of dmps to plot for heatmap",
			"celllines=s@",  "cellines",
			"type=s",	"type of analysis [differential|variability]",
			"niter=i",	"number of iteration",
			"dmps=s@",  "list of dmps id to plot with comet",
			"dmrs=s@",  "list of dmrs id to plot with comet",
			"max=i",    "number of dmps to plot with comet",
			"win=i",     "windows size for comet plot",
			"separator=s", "separator character in files"
			
		), ncol=2, byrow=TRUE))




path <- dirname(strsplit(commandArgs()[4],"=")[[1]][2])
datadir<-paste0(path , "/../data/")
if ( !file.exists(out)) { dir.create(out) }

if(tool=="loader"){
  
  source( paste0(path, "/loader.r")  )
  source( paste0(path, "/utils.r")  )
  source( paste0(path, "/plots.r")  )
  source( paste0(path, "/annot.r")  )
  source( paste0(path, "/pca.r")  )
  source( paste0(path, "/missingValues.r")  )
  
  samples<-suppressWarnings( autoformat(samples) )
  barcode<-suppressWarnings( autoformat(barcode) )
  groups <-suppressWarnings( autoformat(groups) )
  if( filters!="default" ){ filters<-unlist(strsplit(filters,",")) }
  loader()
  
}

if(tool=="batchcorrection"){
  
  source( paste0(path, "/batchcorrection.r")  )
  source( paste0(path, "/utils.r") )
  source( paste0(path, "/pca.r") )
  
  if ( !file.exists(meth)) { 
    stop("Meth file is required. Use 'methylkey -T loader' to create meth file.")
  }
  if (length(variables)==0){
    stop("You must provide at least one variable to protect, usually the sample group ( --variables )")
  }
  print("load data")
  load(meth)
  batch<-suppressWarnings( autoformat(batch) )
  batchcorrection()
  
}

if(tool=="dmps"){
  
  source( paste0(path, "/utils.r")  )
  source( paste0(path, "/plots.r")  )
  source( paste0(path, "/annot.r")  )
  source( paste0(path, "/dmps.r")  )
  
  case=gsub("__cn__","",case) #chomp endline when the class parameter is the last one 
  control=gsub("__cn__","",control) #chomp endline when the class parameter is the last one
  celllines<-unlist(strsplit(celllines,","))
  
  print("load data")
  load(meth)
  dmps()
}


if(tool=="dmrcate"){
  
  source( paste0(path, "/utils.r")  )
  source( paste0(path, "/plots.r")  )
  source( paste0(path, "/annot.r")  )
  source( paste0(path, "/dmrcate.r")  )
  
  print("load data")
  load(meth)
  dmrcate_wrapper()
}


if (tool=="comet"){
  
  source( paste0(path, "/utils.r")  )
  source( paste0(path, "/plots.r")  )
  source( paste0(path, "/annot.r")  )
  source( paste0(path, "/comet_plot.r")  )
  
  print("load data")
  if ( max==0 & length(dmps)==0) { stop("no dmps to plot, set --max or --dmps options") }
  load(meth)
  comet_plot()
  
}

if (tool=="plot_dmrs"){
  
  source( paste0(path, "/utils.r")  )
  source( paste0(path, "/plots.r")  )
  source( paste0(path, "/annot.r")  )
  source( paste0(path, "/dmrs_plot.r")  )
  
  print("load data")
  if ( max==0 & length(dmrs)==0) { stop("no dmrs to plot, set --max or --dmrs options") }
  load(meth)
  plot_dmrs()
  
}

if (tool=="report"){
  
  library(rmarkdown)
  render(paste0(path,"/methylkey.Rmd"), params = list(project=out) )
  
}

