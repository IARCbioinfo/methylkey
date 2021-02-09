message("Welcome to methylkey v2.2")
#0- libraries
message("loading libraries ...")
message(getwd())

source_https <- function(url, ...) {
  # load package
  require(RCurl)
  
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

suppressPackageStartupMessages( library(data.table) )
#suppressPackageStartupMessages( library(knitr) )
suppressPackageStartupMessages( library(wateRmelon) )
suppressPackageStartupMessages( library(minfi) )
suppressPackageStartupMessages( library(GetoptLong) )
suppressPackageStartupMessages( library(EpiSmokEr) )
#IlluminaHumanMethylationEPICmanifest
#IlluminaHumanMethylationEPICanno.ilm10b4.hg19

#path="methylkey/R/"
#source(paste0(path,"sampleSheet.r"))
#source(paste0(path,"missingValues.r"))
#source(paste0(path,"batchcorrection.r"))
#source("methylkey/R/utils.r")

githome="http://git.iarc.lan/EGE/methylkey/raw/master/"
suppressMessages( source_https(paste0(githome,"/R/utils.r")) )
suppressMessages( source_https(paste0(githome,"/R/sampleSheet.r")) )
suppressMessages( source_https(paste0(githome,"/R/missingValues.r")) )
#suppressMessages( source_https(paste0(githome,"/R/batchcorrection.r")) )
source("methylkey/R/batchcorrection.r")
source("methylkey/R/pca.r")

message("Reading parameters ...")
#1 parameters
opt = new.env()
opt$pdata="pdata.csv"
opt$idat="idat"
opt$normalize="Funnorm"
opt$nosva=FALSE
opt$samples=NULL
opt$barcode=NULL
opt$groups="Sentrix_ID"
opt$model="~Smoking_status"
opt$nalimit=0.2
opt$pval=0.02
opt$missing="mean"
opt$separator="\t"
opt$out="minfi_result"
opt$meth=paste0(opt$out,"/betas.rdata")
opt$cell=FALSE
opt$badsamples=""
opt$filter="default"

GetoptLong(
  "pdata=s",      "pdata file",
  "idat=s",       "idat dir",
  "meth=s",       "input object",
  "normalize=s",  "(Funnorm|Illumina|Noob|Quantile|SWAN|SWANoob|None)",
  "sva!",         "batch correction with sva",
  "samples=s",    "sample column name or index",
  "barcode=s",    "barcode column name or index",
  "groups=s",     "groups columns name or index (comma separated list without space)",
  "separator=s",  "file separator",
  "nalimit=f",    "maximum perc. of NA value for each sample",
  "pval=f",       "maximum pvalue for each probe",
  "missing=s",    "missing value method (mean|impute|keep)",
  "filter=s",     "files : list of cpg to remove",
  "badsamples=s", "file : list of sample to remove",
  "cell!",        "houseman",
  "model=s",      "model",       
  "out=s",        "out file name",
  envir=opt
)

if ( !dir.exists(opt$out) ){ dir.create(opt$out) }
groups<-unlist(strsplit(opt$groups,","))
analyse=new.env()

if (opt$filter=="default"){
  opt$filter="http://git.iarc.lan/EGE/methylkey/raw/master/data/Crossreactive_probes_EPIC.csv,http://git.iarc.lan/EGE/methylkey/raw/master/data/SNP_EPIC.csv,http://git.iarc.lan/EGE/methylkey/raw/master/data/SNP_ALL_races_5percent.csv,http://git.iarc.lan/EGE/methylkey/raw/master/data/Sex_EPIC.csv"
}

#################################################
#2- load data
if ( ! file.exists( paste0(opt$out, "/RGset.rda") ) ){
  
  message("loading pdata ...")
  #################################################
  #2.1- load SampleSheet
  pdata<-readSampleSheet(opt$pdata, samples=opt$samples, barcode=opt$barcode, groups=groups, sep=opt$separator)
  if(file.exists(opt$badsamples)){
    pdata<-removeBadSamples(pdata,opt$bad)
  }
  analyse$pdata<-pdata
  analyse$nsamp<-nrow(pdata)

  message("loading idats ...")
  ################################################
  #2.2- read idat
  RGset<-read.metharray.exp(base = opt$idat, targets=data.frame(pdata), force=TRUE)
  betas<-getBeta(RGset)
  #colnames(betas)<-pdata$samples
  analyse$platform<-RGset@annotation["array"]
  message(paste0("reading idat : " , analyse$platform ) )
  analyse$nbprobes1<-nrow(betas)
  message(paste0("reading idat : " , analyse$nbprobes1, " probes" ) )
  analyse$pvalues<-minfi::detectionP(RGset)

  message("predict sex and age ...")
  ################################################
  #2.3- check gender and age
  GMsetEx <- mapToGenome(RGset) 
  estSex <- getSex(GMsetEx)
  pdata$predictedSex <- estSex$predictedSex
  analyse$estSex<-estSex
  #pdata$predicted_BMI_1<-BMIp(betas)
  
  data(coef)
  pdata$predictedAgeHorvath <- agep(betas,method = "horvath")
  data(hannumCoef)
  pdata$predictedAgeHannum <- agep(betas,method = "hannum", coeff=hannumCoef)
  
  #message("Estimate cell count ...")
  ################################################
  #RGset.450k = convertArray(RGset, outType = "lluminaHumanMethylation450k")
  #estimateCellCounts2(RGset, compositeCellType = "Blood",
  #                   processMethod = "auto", probeSelect = "auto",
  #                   cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu"),
  #                   referencePlatform = analyse$platform,
  #                   returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)
  
  message("saving RGset raw ...")
  #####################################################
  #2.4- saving : RGset
  opt_=opt
  save(opt_, analyse, RGset, file=paste0(opt$out, "/RGset.rda") )
}

#########################################################
#3- QC before processing
if ( !file.exists( paste0(opt$out, "/plotQC1.jpg") )){
  
  if (!exists("RGset")){
    message("loading RGset ...")
    load( paste0(opt$out, "/RGset.rda") )
    opt<-appendEnv(opt_,opt)
  }
  
  message("drawing boxplots ...")
  #####################################################
  #3.1- QC before processing : plotQC1
  nsamp=analyse$nsamp
  
  jpeg(paste0(opt$out, "/boxplot_colour1.jpg"), width=800, height=800)
  ylab<-"log2 intensity of both green and red channel"
  par(xaxt='n')
  boxplot(log2(getRed(RGset)+1), col = "red", boxwex = 0.25, at= 1:nsamp - 0.175, ylab=ylab, labels=pdata$samples, cex=0.5)
  boxplot(log2(getGreen(RGset)+1), col = "green", boxwex = 0.25, at= 1:nsamp + 0.175, axis=F , add=T, cex=0.5)
  par(xaxt='s')
  axis(1, at=1:nsamp, labels=pdata$samples, tick=TRUE, las=2, cex.axis=0.8)
  dev.off()
  
  jpeg(paste0(opt$out, "/redchannel.jpg"), width=800, height=800)
  boxplot(log2(getRed(RGset)+1), main="Red channel", col = "red", las=2, cex.axi=0.8, labels=pdata$samples)
  dev.off()
  
  jpeg(paste0(opt$out, "/greenchannel.jpg"), width=800, height=800)
  boxplot(log(getGreen(RGset)+1), main="Green channel", col = "green", las=2, cex.axi=0.8, labels=pdata$samples)
  dev.off()
  
  message("drawing QC plot ...")
  #####################################################
  #3.2- QC before processing : plotQC1
  MSet <- preprocessRaw(RGset)
  MSet <- fixMethOutliers(MSet)
  qc <- getQC(MSet)
  jpeg(paste0(opt$out, "/plotQC1.jpg"), width=800, height=800)
  plotQC(qc)
  dev.off()
}

message("normalization ...")
#########################################################
#4- Normalization
if ( !file.exists( paste0(opt$out, "/", opt$normalize , "/betas.rda") )){
  
  dir.create(paste0(opt$out, "/", opt$normalize), showWarnings = F)
  if (!exists("RGset")){
    message("loading RGset ...")
    load( paste0(opt$out, "/RGset.rda") )
    opt<-appendEnv(opt_,opt)
    estSex<-analyse$estSex
    pdata<-analyse$pdata
    betas<-getBeta(RGset)
    groups<-unlist(strsplit(opt$groups,","))
  }

  message(opt$normalize)
  ######################################################
  #4.1- Normalization
  isna<-which(is.na(betas))
  if(opt$normalize=="Funnorm") {	MSet<-preprocessFunnorm(RGset, sex=estSex$predictedSex)   }
  if(opt$normalize=="Illumina"){	MSet<-preprocessIllumina(RGset)                           }
  if(opt$normalize=="Noob")    {	MSet<-preprocessNoob(RGset)	                              }
  if(opt$normalize=="Quantile"){	MSet<-preprocessQuantile(RGset, sex=estSex$predictedSex)  }
  if(opt$normalize=="SWANoob") {  MSet<-preprocessSWAN(RGset, preprocessNoob(RGset))        }
  if(opt$normalize=="SWAN")    {  MSet<-preprocessSWAN(RGset, preprocessRaw(RGset))         }

  betas=getBeta(MSet)
  betas[ !isna ]<-NA #After normalization NA values are replace by values close to 0. This restore the NA status.
  
  message("Filtering by pvalue ...")
  ######################################################
  #4.2- Filter by pvalue
  betas[ analyse$pvalues[rownames(betas),] > opt$pval ] <- NA
  
  message("get probes to remove from list ...")
  ######################################################
  #4.3- select probes from list
  probes=c()
  filters<-split(",",opt$filter)
  for (file in filters){
    if(file.exists(file)){
      probes <- unique( c(probes, fread(file)[[1]] ))
    }
  }
  filteredFromList<-length(probes)
  message(filteredFromList)
  
  message("get probes to remove with percentage of NA values > nalimit ...")
  ######################################################
  #4.4- Select probes with percentage of NA values > nalimit
  naprobes<-CpGNAexcl( betas,opt$nalimit )
  filteredNAprobes<-length(naprobes)
  write.table(paste(naprobes, collapse="\n"), file=paste0(opt$out, "/removed.txt"), row.names=F, sep="\t")
 
  message("remove selected cpg ...")
  ######################################################
  #4.5- remove selected probes
  probes<-unique( c(probes, naprobes) )
  message( paste0( "total number of removed probes : ", length(probes) ))
  betas<-betas[ ! rownames(betas) %in% probes,]
  
  message("process missing values ...")
  ######################################################
  #4.6- remove remaining missing values
  if ( sum(is.na(betas)) > 0 ) {
    if (opt$missing=="mean"){
      betas<-replaceByMean(betas,pdata[,groups[1]])
    }
    if (opt$missing=="impute"){
      betas<-imputeNA(betas,nalimit)
    }
  }
  
  
  message("predict smoking statut")
  ######################################################
  #4.7- predict smoking statut
  print(pdata$predictedSex)
  pdata$sex<-pdata$predictedSex
  print(pdata$sex)
  colnames(betas)<-rownames(pdata)
  
  analyse$result_SSt <- epismoker(dataset=betas, samplesheet = pdata, method = "SSt")
  pdata$PredictedSmokingStatus <- analyse$result_SSt$PredictedSmokingStatus
  pdata = pdata[,!(names(pdata) %in% "sex")]
  
  ######################################################
  #4.8- pca
  analyse$betas_pca<-makepca( betas, pdata, colnames(pdata), nPC=10 )
  
  ######################################################
  #4.9- Save
  message("saving ...")
  message( paste0(opt$out, "/", opt$normalize ,"/betas.rda") )
  opt_=opt
  save(opt_, analyse, betas, pdata, file=paste0(opt$out, "/", opt$normalize ,"/betas.rda") )
  
}

#########################################################
#5- batch correction
if ( !opt$nosva ){
  
  message("loading betas ...")
  if (!is.matrix("betas")){
    message("loading meth.rda ...")
    load( paste0(opt$out, "/", opt$normalize , "/betas.rda") )
    opt<-appendEnv(opt_,opt)
  }
  
  message("running sva ...")
  ######################################################
  #5.1- batchcorrection sva
  mval<-beta2m(betas)
  mval[!is.finite(mval)]<-min(mval[is.finite(mval)])
  mval<-bc_sva(mval,pdata,opt$model)
  
  ######################################################
  #5.2- pca
  analyse$mval_pca<-makepca( mval, pdata, colnames(pdata), nPC=10 )  
  print(dim(analyse$mval_pca))
    
  ######################################################
  #5.3- Save
  message("saving ...")
  message( paste0(opt$out, "/", opt$normalize ,"/mval.rda") )
  opt_=opt
  save(opt_, analyse, mval, pdata, file=paste0(opt$out, "/", opt$normalize ,"/mval.rda") )
}


message("the end")