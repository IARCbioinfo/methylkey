#0- libraries
message("loading libraries ...")

source_https <- function(url, ...) {
  # load package
  require(RCurl)
  
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

suppressPackageStartupMessages( library(GetoptLong) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(limma) )
suppressPackageStartupMessages( library(qqman) )
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(wateRmelon) )

path="methylkey/R/"
#source(paste0(path,"/utils.r"))
#source(paste0(path,"/regression.r"))

githome="http://git.iarc.lan/EGE/methylkey/raw/master/"
suppressMessages( source_https(paste0(githome,"/R/utils.r")) )
suppressMessages( source_https(paste0(githome,"/R/regression.r")) )

message("Reading parameters ...")
#1 parameters
opt=new.env()
opt$meth="minfi_result/Funnorm/betas.rda"
opt$model="~group"
opt$case="case"
opt$control="control"
opt$fdr=0.05
opt$pval=0.1
opt$reg="ls"
opt$niter=25
opt$ncore=4
opt$genome="hg19"
opt$dmrcate=FALSE
opt$annot=FALSE
opt$annotatr=NULL
opt$fullmanifest=FALSE

GetoptLong(
  "meth=s",       "input object",
  "model=s",      "model",
  "case=s",       "case",
  "control=s",    "control",
  "fdr=f",        "fdr threshold for DMPs, default=0.05",
  "pval=f",       "pvalue threshold for DMRs, default=0.1",
  "reg=s",        "method of regression, default=ls",
  "niter=i",      "number of iteration, default=25",
  "ncore=i",      "number of core, default=4",
  "fullmanifest!","annotate probes with all columns of the manifest",
  "dmrcate!",     "calculate DMRs",
  "annot!",       "annot DMRs",
  "genome=s",     "ref genome, default hg19",
  "annotatr=s",   "use prebuild annotation, default NULL",
  envir=opt
)
  
message("loading data ...")
#################################################
#1.1- load data
load(opt$meth)
message(analyse$platform)
if(file.exists(gsub("betas","mval",opt$meth))){ load(gsub("betas","mval",opt$meth)) }
message(analyse$platform)
if(!exists("betas")){ stop("betas not found ! check your input file") }
message(analyse$platform)
if(!exists("mval")){
  mval<-beta2m(betas) 
  mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
}

message("debug")
message(analyse$platform)
message(opt$genome)
message(paste0("number of probes in mvalues : ",nrow(mval)))

opt<-appendEnv(opt_,opt)                 # update options with new values
message(analyse$platform)
message(opt$genome)
model=gsub("~","",opt$model)
rpath<-paste(dirname(opt$meth), model , sep="/") #result path
dir.create(rpath,recursive=TRUE, showWarnings = FALSE)
print(rpath)
log=new.env()
log$comp=paste(opt$case,"vs",opt$control,sep="-")

if (!file.exists( paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmps.rda") )){ 
  
  #1.2- remove cpgs when sum of its betas is 0
  #it can produce a BUG when it is the case
  iszero<-which(rowSums(mval) == 0 )
  if(length(iszero) > 0){
    mval<-mval[-iszero,]
  }
  
  message("building design ...")
  #################################################
  #2.1- levels
  case<-make.names(opt$case)
  control<-make.names(opt$control)
  print(case)
  print(control)
  variables=unlist(strsplit(gsub("~","",opt$model),"\\+"))
  print(variables[1])
  if (sum( ! variables[1] %in% colnames(pdata) ) ){ stop("your model do not match samplesheet columns names") }
  print( as.factor(pdata[,variables[1]]) )
  #lvls<-levels(as.factor(pdata[,variables[1]]))
  lvls<-unique( unlist( pdata[,variables[1]] ) )
  print(lvls)
  #levels<-paste0(levels(as.factor(pdata[,variables[1]])),collapse=",")
  #save.image("DEBUG.rda")
  if(! case %in% lvls ){ stop( paste0( "Error : ", case , " not in ",  lvls ) )}
  if(! control %in% lvls ){ stop( paste0( "Error : ", control , " not in (", paste(lvls),")" )  )}
  pdata[,variables[1]] <- relevel(as.factor(unlist(pdata[,variables[1]])), case)
  pdata[,variables[1]] <- relevel(as.factor(unlist(pdata[,variables[1]])), control)
  message(opt$model)
  message(log$comp)
  message(opt$model)
  #2.2- design
  design<-model.matrix(as.formula(opt$model),data=pdata)
  colnames(design)<-make.names(colnames(design))
  cmtx <- makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )
  
  message("Running regression ...")
  #################################################
  #2.1- Run regression analysis for dmps
  save.image(file="DEBUG.rda")
  regression<-m_regression(mval,pdata,variables,method=opt$reg,niter=opt$niter,ncore=opt$ncore)
  log$lambda=regression$lambda
  
  message(paste0("number of probes in regression result", nrow(regression$table)))

  message("get delta betas ...")
  #################################################
  #3.1- Get DeltaBetas
  message("debug")
  regression$table$deltaBetas<-getDeltaBetas(betas[rownames(regression$table),],pdata[,variables[1]],case,control)*100
  print(dim(regression$table))
  
  message("QQplot ...")
  #################################################
  #4.1- Draw qqplot
  jpeg( paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, "_qqplot.jpg") )
  try(
    qq(regression$pvals,main=paste("QQ plot: ","lambda=", regression$lambda,sep=""))
  ,silent=TRUE)
  dev.off()
  
  message("Downloading annotation ...")
  #################################################
  message(analyse$platform)
  message(opt$genome)
  if( grepl("27k",analyse$platform)  & opt$genome=="hg19"){ manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM27/HM27.hg19.manifest.tsv.gz") }
  if( grepl("450k",analyse$platform) & opt$genome=="hg19"){ manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM450/HM450.hg19.manifest.tsv.gz") }
  if( grepl("EPIC",analyse$platform) & opt$genome=="hg19") { manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/EPIC/EPIC.hg19.manifest.tsv.gz") }
  if( grepl("27k",analyse$platform)  & opt$genome=="hg38"){ manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM27/HM27.hg38.manifest.tsv.gz") }
  if( grepl("450k",analyse$platform) & opt$genome=="hg38"){ manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz") }
  if( grepl("EPIC",analyse$platform) & opt$genome=="hg38") { manifest<-fread("https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.tsv.gz") }
  if( grepl("MM285",analyse$platform) ) {
    library(dplyr)
	source_https("http://git.iarc.lan/EGE/methylkey/raw/master/R/annotation.r")
    #manifest<-fread("https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/Infinium%20Mouse%20Methylation%20v1.0%20A1%20GS%20Manifest%20File.csv",skip=7) 
    manifest<-fread("https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv")
    #manifest<-manifest[,c("chrom","chromStart","chromEnd","chromLength","chromStrand","name","Gene","Transcript","Feature","Source","CpG_Island_chromStart","CpG_Island_chromEnd","CpG_Island_chromEnd","CpG_Island_cpgNum","CpG_Island_gcNum","CpG_Island_perCpg","CpG_Island_perGc","CpG_Island_obsExp","S_Shore","S_Shelf","MFG_Change_Flagged")]
    colnames(manifest)<-c("probeID","Gene","Transcript","chr","start","end","length","strand","Feature","Source")
	#manifest<-makeGRangesFromDataFrame(manifest[!grepl("rs",manifest$probeID) & manifest$chr>0,], keep.extra.columns=TRUE)
        #seqlevelsStyle(manifest)<-"UCSC" 
        #annotated_manifest<-mk_annotatr(manifest, genome="mm10")
	#foo<-unique(as.data.frame(annotated_manifest)[,c(1:3,5,6,30)])
	#foo<-foo %>% arrange(probeID,annot.symbol)
	#p <- function(v) { Reduce(f=paste, x = v) }
	#manifest <- foo %>% group_by(seqnames,start,end,strand,probeID) %>% summarize( genes=p(annot.symbol) ) %>% ungroup() %>% as.data.table()
        #manifest$genes <- gsub(" *NA","",manifest$genes)
  }

  message("annotating ...")
  save.image("debug.rda")
  #################################################
  if(opt$fullmanifest){
    table<-merge(manifest[,1:ncol(manifest)],data.table(probeID=rownames(regression$table), regression$table), by="probeID")
  } else if ( grepl("MM285",analyse$platform) ) {
    foo<-data.table(probeID=rownames(regression$table), regression$table)
    table<-merge(manifest, data.table(probeID=rownames(regression$table), regression$table), by="probeID")
  } else {
    print(dim(regression$table))
    table<-merge(manifest[,c(1:5,20)],data.table(probeID=rownames(regression$table), regression$table), by="probeID")
    print(dim(regression$table))
  }
  colnames(table)[2:5]<-c("chr","start","end","strand")
  
  message("saving ...")
  #################################################
  #5.1- Save
  opt_=opt
  write.table(table,sep="\t",row.names=FALSE,quote=F,file=paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmps.tab"))
  save(opt_, log, table, file=paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmps.rda") )
}

if (opt$dmrcate){
  
  suppressPackageStartupMessages( library(DMRcate) )
  suppressPackageStartupMessages( library(DMRcatedata) )
  
  message("DMR analysis ...")
  ################################################
  if(!exists("regression")){
    message("Try loading DMPs...")
    load(paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmps.rda"))
    opt<-appendEnv(opt_,opt)                 # update options with new values
    #rpath<-paste(opt$out, opt$normalize, model, sep="/") #result path
  }
  
  message("DMRcate ...")
  ###############################################
  table<-makeGRangesFromDataFrame(table, keep.extra.columns=TRUE)
  annotated <- data.frame(chr=seqnames(table), start=start(table), end=end(table), strand=strand(table),
                          stat=table$t, diff= table$deltaBetas, ind.fdr=table$adj.P.Val, is.sig=(table$adj.P.Val<opt$fdr) )
  annotated<-makeGRangesFromDataFrame(annotated, keep.extra.columns=TRUE)
  names(annotated)<-table$probeID
  foo <- new("CpGannotated", ranges=sort(annotated))
  message(length(foo@ranges))
  message(length(foo@ranges$stat))
  if( sum(is.na(foo@ranges$diff)) ){foo@ranges$diff[ which(is.na(foo@ranges$diff)) ] <- 0 } # fix issue when there is missing values in betas.

  dmrcoutput<-dmrcate(foo,C=2, pcutoff=opt$pval)
  table <- extractRanges(dmrcoutput, genome = opt$genome)
  
  message("Retrieve overlaping probes ...")
  ###############################################
  #retrieve overlaping sites (probes)
  overlap <- findOverlaps(table,annotated)
  overlap <- data.table( cbind( query=queryHits(overlap), subject=subjectHits(overlap) ) )
  table$overlapping.sites<-overlap[ , paste( names(annotated)[subject], collapse="," ) , by=query ]$V1

  message("saving ...")
  #################################################
  #5.1- Save
  opt_=opt
  write.table(as.data.frame(table),sep="\t",row.names=FALSE,quote=F,file=paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmrs.tab"))
  save(log, opt_, table, file=paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmrs.rda") )
  
}


if (opt$annot){
  
  source_https("http://git.iarc.lan/EGE/methylkey/raw/master/R/annotation.r")
  
  if(!is.data.table(table)){
    message("Try loading DMRs...")
    load(paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".dmrs.rda"))
    opt<-appendEnv(opt_,opt)                 # update options with new values
    #rpath<-paste(opt$out, opt$normalize, model, sep="/") #result path
  }
  
  message("Annotating DMRs ...")
  #################################################
  #6.1- Annot dmrs
  annot=NULL
  if(!opt$genome %in% c("hg19","hg38","mm10")){ message(paste0(opt$annot," is not a valid genome")) }
  if (!is.null(opt$annotatr)){load(opt$annotatr)}
  
  seqlevelsStyle(table)<-"UCSC" 
  table<-mk_annotation(regions=table,annot=annot,genome=opt$genome)
  message(paste0("table:",length(table)))
  message("saving ...")
  #################################################
  #6.2- Save
  opt_=opt
  print(paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".annotated.dmrs.rda"))
  save(table, log, opt_, file=paste0(rpath, "/", model, "_", log$comp, "_", opt$reg, ".annotated.dmrs.rda") )
  
}

message("The end ...")
