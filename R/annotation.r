#####################################################################################
#
# Title  : annotation.r
# Author : CahaisV@iarc.fr
# Package: methylkey
# Date   : 02/07/2020
#
#####################################################################################

library(data.table)
library(GenomicRanges)
library(annotatr)

#example : annotation of dmrs
#load("dmrs.rda")
#load("../methylkey/annotatr/annot.hg19.rda")
#annotated_dmrs<-mk_annotation(table,annot)

#example of annotation from a gtf file
#gr<-gtf2GenomicRanges("hg19.ensGene.gtf.gz")
#gr<-mk_annotation(gr,genome="hg19")

#example : create GenomeInfo object from UCSC
#hg38<-GenomeInfoDb::fetchExtendedChromInfoFromUCSC("hg38")
#genInfo<-GenomeInfoDb::Seqinfo(hg38$UCSC_seqlevel, hg38$UCSC_seqlength, hg38$circular, "hg38")
#genInfo<-keepStandardChromosomes(genInfo)

#example : annotate 850k manifest
#manifest<-fread("illumina/EPIC.hg19.manifest.tsv.gz")
#colnames(manifest)[1:4]<-c("chr","start","end","strand")
#manifest<-makeGRangesFromDataFrame(manifest[!is.na(chr),], keep.extra.columns=TRUE)
#annotated_manifest<-mk_annotatr(manifest, genome="hg19")
#annotated_manifest<-as.data.table(annotated_manifest)
#save(annotated_manifest,file="illumina/EPIC.hg19.annotatr.rda")

#######################
#Select tool for annotation
mk_annotation<-function(regions,genome="hg19",annot=NULL,tool="annotatr"){
  
  if (tool=="annotatr") { cpg_annotated=mk_annotatr(regions=regions,annot=annot,genome=genome) }
  return(cpg_annotated)
}


#######################
# build annotation (for annotatr)
buildAnnot<-function(genome="hg19"){
  
  listOfAnnotations<-builtin_annotations()
  annots<-listOfAnnotations[grep(genome, listOfAnnotations)]
  annotations <- build_annotations(genome = genome, annotations = annots)
  return(annotations)
}

#######################
# make annotation using annotatr package
mk_annotatr<-function(regions,annot=NULL,genome="hg19"){
  
  if(is.null(annot)){ annot<-buildAnnot(genome=genome) }
  
  foo<-elementMetadata(annot)
  foo<-as.data.table(foo)
  foo[, c("genome","feature","value"):=tstrsplit(type,"_") ]
  elementMetadata(annot)<-foo
  
  #n<-nrow(foo[ feature=="cpg" & value=="shores", ])
  #shores<-rep(c("northShores","southShores"), n/2)
  #foo[ feature=="cpg" & value=="shores",]$value="shores"
  
  #n<-nrow(foo[ feature=="cpg" & value=="shelves", ])
  #shelves<-rep(c("northShelves","southShelves"), n/2)
  #foo[ feature=="cpg" & value=="shelves",]$value="shelves"
  cpg_annotated <- annotate_regions(regions, annotations=annot, minoverlap = 1L, ignore.strand = FALSE)
  
  return(cpg_annotated)
  
}

#######################
#Read a gtf file and make a Grange object from it.
gtf2GenomicRanges<-function(file){
  
  gtf<-fread(file)
  colnames(gtf)<-c("Seqname","source","feature","start","end","score","strand","frame","attributes")
  gr<-makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE)
  #TODO : set the  seqinfo in gr object
  #seqinfo<-fetchExtendedChromInfoFromUCSC("hg19")
  return(gr)
}

