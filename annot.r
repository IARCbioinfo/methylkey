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
	library(pkg, character.only = TRUE)
	txdb=eval(parse(text = pkg))
	return(txdb)
}


#######################
# build annotation (annotatr)
buildAnnot<-function(genome="hg19"){

	require(annotatr)
	listOfAnnotations<-builtin_annotations()
	annots<-listOfAnnotations[grep(genome, listOfAnnotations)]
	annotations <- build_annotations(genome = genome, annotations = annots)
	return(annotations)
}

#######################
# read regions
readRegions<-function(regions, genome){

	require(annotatr)
	if (is.character(regions)){
		cpg_regions <- read_regions(regions, genome=genome, format="BED")
	} else { 
		cpg_regions <- regions 
	}

}


#######################
# annotation (annotatr)
mk_annotate<-function(regions,genome="hg19"){
	
	require(annotatr)
	annotations=NULL
	if (genome=="hg19"){ load(paste0(path,"/data/annotatr.hg19.rda")) }
	else if (genome=="hg38"){ load(paste0(path,"/data/annotatr.hg38.rda")) }
	else { annotations<-buildAnnot(genome=genome) }

	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shores")])
	shores<-rep(c(paste0(genome,"_cpg_southShores"), paste0(genome,"_cpg_northShores")), (n+1)/2)
	suppressWarnings(annotations[annotations$type==paste0(genome,"_cpg_shores")]$type<-shores[1:n])

	n<-length(annotations[annotations$type==paste0(genome,"_cpg_shelves")])
	shelves<-rep(c(paste0(genome,"_cpg_southShelves"), paste0(genome,"_cpg_northShelves")), (n+1)/2)
	suppressWarnings(annotations[annotations$type==paste0(genome,"_cpg_shelves")]$type<-shelves[1:n])

	cpg_regions <- readRegions(regions, genome)
	cpg_annotated <- annotate_regions(cpg_regions, annotations=annotations, minoverlap = 0L)
	return(cpg_annotated)
}


