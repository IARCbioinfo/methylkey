#####################################################################################
#
# Title  : annot.r
# Author : CahaisV@iarc.fr
# Package: methylkey
# Date   : 12/07/2018
#
#####################################################################################

#######################
#get genome annotation
getTxdb<-function(genome)
{
	if      (grepl("hg",genome)) { species="Hsapiens"  }
	else if (grepl("mm",genome)) { species="Mmusculus" }
	else { return ("no reference") }
	
	suppressMessages(pkg<-paste("TxDb", species, "UCSC", genome ,"knownGene", sep="."))
	print(pkg)
	suppressPackageStartupMessages( library(pkg, character.only = TRUE) )
	txdb=eval(parse(text = pkg))
	return(txdb)
}


#######################
# build annotation (annotatr)
buildAnnot<-function(genome="hg19"){

	suppressPackageStartupMessages(require(annotatr))
	listOfAnnotations<-builtin_annotations()
	annots<-listOfAnnotations[grep(genome, listOfAnnotations)]
	annotations <- build_annotations(genome = genome, annotations = annots)
	return(annotations)
}

#######################
# read regions
readRegions<-function(regions, genome){

	suppressPackageStartupMessages(require(annotatr))
	if (is.character(regions)){
		suppressMessages( cpg_regions <- read_regions(regions, genome=genome, format="BED") )
	} else { 
		cpg_regions <- regions 
	}

}


#######################
# annotation (annotatr)
getRegionsAsGenomicRanges<-function(regions=NULL, platform="array"){
  
  if (is.null(regions)){stop("\r getRegionsAsGenomicRanges : regions cannot be NULL !")}
  if (typeof(regions)!="S4"){
    if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
    if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
    if(regions == "") { stop( paste0("\r annoTatble : Oups ! no regions to annotate : ", platform) )}
    #return(makeGRangesFromDataFrame(regions))
    return(regions)
  }
  return(regions)
}


mk_annotate<-function(regions,genome="hg19"){
	
	suppressPackageStartupMessages(require(annotatr))
	annotations=NULL
	if (genome=="hg19"){ load(paste0(datadir,"/annotatr.hg19.rda")) }
	else if (genome=="hg38"){ load(paste0(datadir,"/annotatr.hg38.rda")) }
	else { annotations<-buildAnnot(genome=genome) }
		
	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shores")])
	shores<-rep(c(paste0(genome,"_cpg_southShores"), paste0(genome,"_cpg_northShores")), (n+1)/2)
	suppressWarnings( annotations[annotations$type==paste0(genome,"_cpg_shores")]$type<-shores[1:n] )

	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shelves")])
	shelves<-rep(c(paste0(genome,"_cpg_southShelves"), paste0(genome,"_cpg_northShelves")), (n+1)/2)
	suppressWarnings( annotations[annotations$type==paste0(genome,"_cpg_shelves")]$type<-shelves[1:n] )

	cpg_regions <- readRegions(regions, genome)
	suppressMessages( cpg_annotated <- annotate_regions(cpg_regions, annotations=annotations, minoverlap = 0L) )
	return(cpg_annotated)
}


annoTable<-function(regions=NULL, toptable=NULL, genome="hg19")
{
  if (is.null(regions)){stop("\r annoTable : regions cannot be NULL !")}
  if (is.null(toptable)){stop("\r annoTable : toptable cannot be NULL !")}
  if (ncol(toptable)!=7){stop("\r annoTable : toptable should have 7 columns !")}
  
  #if (typeof(regions)!="S4"){
  #  if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
  #  if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
  #  if(regions == "") { stop( "\r annoTatble : Oups ! no region to annotate.") }
  #}
  
  #complete annotation
  cpg_annotated<-mk_annotate(regions,genome=genome)
  dm_annotated<-cpg_annotated[ cpg_annotated$name %in% rownames(toptable)]
  dm<-merge(toptable, dm_annotated, by.x="row.names", by.y="name")
  colnames(dm)[1]<-"cpg"
  colnames(dm)[9]<-"chr"
  #second column should be logFC or Estimate
  dm$DM_Status[ as.numeric(dm[,2]) > 0 ] = "hyper"
  dm$DM_Status[ as.numeric(dm[,2]) < 0 ] = "hypo"
  
  #toptable
  dm_annotated<-makeGRangesFromDataFrame(dm, keep.extra.columns=TRUE)
  dm<-dm[order(dm$annot.symbol),]
  toptable<-dm[!duplicated(dm$cpg),c(9:11,1,13,2:8,24)] 
  toptable<-toptable[order(toptable$P.Value),]
  
  return(list(full=dm,toptable=toptable,gr=dm_annotated,cpg_annotated=cpg_annotated))
}
