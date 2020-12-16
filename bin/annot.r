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
	
	require(annotatr)
	annotations=NULL
	if (genome=="hg19"){ load(paste0(datadir,"/annotatr.hg19.rda")) }
	else if (genome=="hg38"){ load(paste0(datadir,"/annotatr.hg38.rda")) }
	else { annotations<-buildAnnot(genome=genome) }
		
	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shores")])
	shores<-rep(c(paste0(genome,"_cpg_southShores"), paste0(genome,"_cpg_northShores")), (n+1)/2)
	annotations[annotations$type==paste0(genome,"_cpg_shores")]$type<-shores[1:n]

	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shelves")])
	shelves<-rep(c(paste0(genome,"_cpg_southShelves"), paste0(genome,"_cpg_northShelves")), (n+1)/2)
	annotations[annotations$type==paste0(genome,"_cpg_shelves")]$type<-shelves[1:n]

	cpg_regions <- readRegions(regions, genome)
	cpg_annotated <- annotate_regions(cpg_regions, annotations=annotations, minoverlap = 0L)
	return(cpg_annotated)
}



# annoTable
# @Desc   : Annotate a Genomic Ranges and return a dataTable
# 			Each line of the regions object is duplicated for each annotation
# @input  : region=GRanges
#			genome=c("hg19","hg38","mm10","rn6")
# @output : DataTable
annoTable<-function(regions, genome="hg19"){

	require(data.table)

	#Annotate
	annotation<-mk_annotate(regions,genome)
	annotation<-as.data.table(annotation)
	return(annotation)
}

#getOfficilaSymbols
# @inputs  : geneList as vector
#		     genome=c("hg19","hg38","mm10","rn6")
# @outputs : geneList as vector
getOfficialSymbols<-function(geneList, genome){

	require(geneSynonym)

	#replace Synonyme gene Symbols by official names
	tax=9606 # homo sapiens
	if (genome == "mm9" || "mm10") { tax=10090 } # Mus musculus
	if (genome == "rn6")           { tax=10116 } # Rattus Norvegicus

	geneList<-unique(unlist(geneList), na.rm=TRUE)
	geneList<-geneSynonym( symbollist, tax)
	return(geneList)
}




##annoTable<-function(regions=NULL, toptable=NULL, genome="hg19")
##{
##  if (is.null(regions)){stop("\r annoTable : regions cannot be NULL !")}
 ## if (is.null(toptable)){stop("\r annoTable : toptable cannot be NULL !")}
 ## if (ncol(toptable)!=7){stop("\r annoTable : toptable should have 7 columns !")}
  
  #if (typeof(regions)!="S4"){
  #  if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
  #  if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
  #  if(regions == "") { stop( "\r annoTatble : Oups ! no region to annotate.") }
  #}
  
  #complete annotation
 ## reg_annotated<-mk_annotate(regions,genome=genome)
 ## dm_annotated<-reg_annotated[ reg_annotated$name %in% rownames(toptable)]
 ## dm<-merge(toptable, dm_annotated, by.x="row.names", by.y="name")
 ## colnames(dm)[1]<-"cpg"
 ## colnames(dm)[9]<-"chr"
  #second column should be logFC or Estimate
 ##dm$DM_Status[ as.numeric(dm[,2]) > 0 ] = "hyper"
 ## dm$DM_Status[ as.numeric(dm[,2]) < 0 ] = "hypo"
  
  #toptable
 ## dm_annotated<-makeGRangesFromDataFrame(dm, keep.extra.columns=TRUE)
 ## dm<-dm[order(dm$annot.symbol),]
 ## toptable<-dm[!duplicated(dm$cpg),c(9:11,1,13,2:8,23,24)] 
  ##toptable<-toptable[order(toptable$P.Value),]
  
  ##return(list(full=dm,toptable=toptable,gr=dm_annotated,reg_annotated=reg_annotated))
##}


annotateSites<-function(table,manifest){

	colnames(manifest)<-c("chr","start","end","probe","V5","strand")
	return( merge( manifest[,c(1,2,3,6,4)], table, by.x="probe", by.y="row.names") )
}