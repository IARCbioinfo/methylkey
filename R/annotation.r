#' buildAnnot
#' 
#' @param genome genome
#' 
#' @return annotation file
#' 
#' @export
#' 
buildAnnot<-function(genome="hg19"){
  
  listOfAnnotations<-builtin_annotations()
  annots<-listOfAnnotations[grep(genome, listOfAnnotations)]
  annotations <- build_annotations(genome = genome, annotations = annots)
  return(annotations)
}

#' make annotation with annotatr
#' 
#' @param regions GRanges
#' @param genome genome
#' @param annot annotation database
#' 
#' @return annotation file
#' 
#' @export
#' 
mk_annotatr<-function(regions,annot=NULL,genome="hg19"){
  
  if(is.null(annot)){ annot<-buildAnnot(genome=genome) }
  
  #foo<-elementMetadata(annot)
  #foo<-as.data.table(foo)
  #foo[, c("genome","feature","value"):=tstrsplit(type,"_") ]
  
  elementMetadata(annot)<-elementMetadata(annot) %>% as.data.frame() %>% tidyr::separate(type,into=c("genome","feature","value"))
  
  cpg_annotated <- annotate_regions(regions, annotations=annot, minoverlap = 1L, ignore.strand = FALSE)
  
  return(cpg_annotated)
}


#' gtf2GenomicRanges
#' 
#' transform gtf file to genomic ranges
#' 
#' @param file gtf file
#' 
#' @return genomic ranges
#' 
#' @export
#' 
gtf2GenomicRanges<-function(file){
  
  gtf<-fread(file)
  colnames(gtf)<-c("Seqname","source","feature","start","end","score","strand","frame","attributes")
  gr<-makeGRangesFromDataFrame(gtf, keep.extra.columns=TRUE)
  #TODO : set the  seqinfo in gr object
  #seqinfo<-fetchExtendedChromInfoFromUCSC("hg19")
  return(gr)
}



#' getInfiniumAnnotation
#' 
#' retrieve infinium annotations
#' 
#' @param plateform plateform
#' @param genome genome
#' 
#' @return data.frame
#' 
#' @export
#' 
getInfiniumAnnotation<-function(plateform, genome="hg19"){

  webAddress="https://zwdzwd.s3.amazonaws.com/InfiniumAnnotation/20180909/"
  webAddress2="https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/"
  
  if( grepl("27k",plateform)  & genome=="hg19") { annot<-read_tsv(paste0(webAddress,"HM27/HM27.hg19.manifest.tsv.gz")) }
  if( grepl("450k",plateform) & genome=="hg19") { annot<-read_tsv(paste0(webAddress,"HM450/HM450.hg19.manifest.tsv.gz")) }
  if( grepl("EPIC",plateform) & genome=="hg19") { annot<-read_tsv(paste0(webAddress,"EPIC/EPIC.hg19.manifest.tsv.gz")) }
  if( grepl("27k",plateform)  & genome=="hg38") { annot<-read_tsv(paste0(webAddress,"HM27/HM27.hg38.manifest.tsv.gz")) }
  if( grepl("450k",plateform) & genome=="hg38") { annot<-read_tsv(paste0(webAddress,"HM450/HM450.hg38.manifest.tsv.gz")) }
  if( grepl("EPIC",plateform) & genome=="hg38") { annot<-read_tsv(paste0(webAddress,"EPIC/EPIC.hg38.manifest.tsv.gz")) }
  if( grepl("MM",plateform)   & genome=="mm10") { 
    annot<-read_csv(paste0(webAddress2,"MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv"))
  }
  
  return(annot)
}

#' getAnnotedManifest
#' 
#' retrieve infinium annotations
#' 
#' @param plateform plateform
#' 
#' @return annotated manifest
#' 
#' @export
#' 
getAnnotedManifest<-function(plateform){
  
  manifest=NULL
  
  if(plateform=="IlluminaHumanMethylation450k"){
    manifest<-readr::read_csv("~/git/methylkey/v1.0/data/HumanMethylation450_15017482_v1-2.csv", skip=7) %>% 
      dplyr::rename( Gene = UCSC_RefGene_Name ) %>% 
      dplyr::rename( chr = CHR ) %>%
      mutate( chr = paste0("chr",chr) ) %>%
      mutate( start = MAPINFO ) %>%
      mutate( end = MAPINFO ) %>%
      mutate( Strand = ifelse(Strand=="F","+","-") )
  }
  
  if(plateform=="IlluminaHumanMethylationEPIC"){
    manifest<-readr::read_csv("~/git/methylkey/v1.0/data/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7) %>% 
      dplyr::rename( probeID = Name ) %>% 
      dplyr::rename( Gene = UCSC_RefGene_Name ) %>% 
      dplyr::rename( chr = CHR ) %>%
      mutate( chr = paste0("chr",chr) ) %>%
      mutate( start = MAPINFO ) %>%
      mutate( end = MAPINFO ) %>%
      mutate( Strand = ifelse(Strand=="F","+","-") )
  }
    
  if(plateform=="IlluminaMouseMethylation285k"){
    manifest<-readr::read_csv("~/git/methylkey/v1.0/data/MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv") %>% 
      dplyr::rename( chr = chrom ) %>%
      dplyr::rename( start = chromStart ) %>%
      dplyr::rename( end = chromEnd ) %>%
      dplyr::rename( Strand = chromStrand )
  }
  
  return(manifest)
}

#' initEnrichR
#' 
#' initiate EnrichR database connexion
#' 
#' @param genome genome
#' 
#' @return dbs
#' 
#' @export
#'
initEnrichR<-function(genome="hg19"){
  
  dbs=NULL
  enrichR::setEnrichrSite("Enrichr")
  websiteLive <- TRUE
  dbs <- enrichR::listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  
  if (websiteLive & genome=="mm10") {
    dbs<-c( dbs$libraryName[grepl("Mouse",dbs$libraryName)], 
              "GO_Biological_Process_2021",
              "GO_Molecular_Function_2021",
              "GO_Cellular_Component_2021",
              "Descartes_Cell_Types_and_Tissue_2021",
              "MGI_Mammalian_Phenotype_Level_4_2021",
              "CellMarker_Augmented_2021")
  }
  
  if (websiteLive & genome %in% c("hg19","hg38")) {
    dbs<-c( "Human_Gene_Atlas", "TRANSFAC_and_JASPAR_PWMs", "KEGG_2021_Human", "GO_Biological_Process_2021",
            "GO_Molecular_Function_2021","GO_Cellular_Component_2021","WikiPathway_2021_Human")
  }
  
  return(dbs)
}