#' violin_plot1
#' 
#' Display betas in a violin plot by group
#' 
#' @param betas matrix of betas
#' @param group variable group
#' @param numPositions number of random CpG used
#' 
#' @return plot
#' 
#' @export
#' 
violin_plot1<-function(betas, group, numPositions=1000){
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-cbind(group=as.character(group),data.table(t(cpg)) )
  cpg<-melt( cpg, id.vars="group" )
  colnames(cpg)<-c("group","probeID","betas")
  
  p1<-ggplot(cpg, aes(x=group, y=betas, fill=group)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("")
  
  return(p1)
}

#' violin_plot2
#' 
#' Display betas in a violin plot by group
#' 
#' @param betas matrix of betas
#' @param annot annotation table
#' @param feature annotation feature to use.
#' @param numPositions number of random CpG used
#' 
#' @return plot
#' 
#' @export
#' 
violin_plot2<-function(betas, annot, feature, numPositions=1000){
  
  require(ggplot2)
  require(data.table)
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-melt(cpg)
  colnames(cpg)<-c("probeID","Basename","betas")
  cpg<-merge(annot,cpg,by="probeID")
  
  #violin plots
  p1<-ggplot(cpg, aes_string(x=feature, y="betas", fill=feature)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("")
  
  return(p1)
}

#' violin_plot3
#' 
#' Display betas in a violin plot by group
#' 
#' @param betas matrix of betas
#' @param group variable group
#' @param annot annotation table
#' @param feature annotation feature to use.
#' @param numPositions number of random CpG used
#' 
#' @return plot
#' 
#' @export
#' 
violin_plot3<-function(betas, annot, group, feature, numPositions=1000){
  
  require(ggplot2)
  require(data.table)
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-cbind(group=as.character(group),data.table(t(cpg)) )
  cpg<-melt( cpg, id.vars="group" )
  colnames(cpg)<-c("group","probeID","betas")
  cpg<-merge(annot,cpg,by="probeID")
  
  #violin plots
  p1<-ggplot(cpg, aes_string(x=feature, y="betas", fill=feature)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("") +
    facet_wrap(~group)
  
  return(p1)
}



#' circosplot
#' 
#' Display betas in a violin plot by group
#' 
#' @param ranges GRanges object
#' @param genome (hg38,hg19,mm10)
#' 
#' @return plot
#' 
#' @export
#' 
circosplot<-function(ranges, genome){

	suppressPackageStartupMessages( require(ggbio) )

	if (genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38) ; species=BSgenome.Hsapiens.UCSC.hg38}
	if (genome=="hg19"){ require(BSgenome.Hsapiens.UCSC.hg19) ; species=BSgenome.Hsapiens.UCSC.hg19}
  if (genome=="hg19"){ require(BSgenome.Mmusculus.UCSC.mm10) ; BSgenome.Mmusculus.UCSC.mm10}

	chr.len = seqlengths(species)
	chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
	myIdeo <- GRanges(seqnames = names(chr.len), ranges = IRanges(start = 1, chr.len))
	seqlengths(myIdeo)<-myIdeo@ranges@width

	#seqlevels(ranges) = names(chr.len)
	ranges<-keepSeqlevels(ranges, names(chr.len), pruning.mode="coarse")
	seqlengths(ranges)<-myIdeo@ranges@width
	
	g.po <-ranges[ranges$meth.diff < 0 ]
	g.per<-ranges[ranges$meth.diff > 0 ]

	p<- ggbio() + circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
	
	if (length(g.po) > 0){
		values(g.po)$id = "hypo"
		p <- p + circle(g.po, geom = "point", size = 1, aes(x = midpoint, y = "meth.diff", color = id), radius = 19, trackWidth = 20) + scale_colour_manual(values = c("magenta", "green")) 
	}
	if (length(g.per) > 0){	
		values(g.per)$id = "hyper"
		p <- p + circle(g.per, geom = "point", size = 1, aes(x = midpoint, y = "meth.diff", color = id), radius = 41, trackWidth = 20)
	}
	# bug in ggbio waiting for fix
	#p<- p + circle(ranges, geom = "text", aes(label = seqnames), vjust = 0, radius = 55, trackWidth = 7)
	return(p)

}




#' barplots_
#' 
#' Display repartition of DMPs/DMRs annotation
#' 
#' @param plateform plateforme name
#' @param tab cpg table
#' @param man manifest
#' 
#' @return plot
#' 
#' @export
#' 
barplots_<-function(plateform="IlluminaHumanMethylation450k",tab=NULL,man=NULL){
  
  cat(paste0("\nLocation annotation for ",plateform," probes, hyper/hypo methylated probes and significantly hyper/hypo methylated probes\n"))
  if( plateform=="IlluminaHumanMethylation450k") { barplots_450k(tab,man) }
  if( plateform=="IlluminaMouseMethylation285k") { barplots_MM280(tab,man) }
  if( plateform=="IlluminaHumanMethylationEPIC") { barplots_EPIC(tab,man) }
                   
}

#' barplots_MM280
#' 
#' Display repartition of DMPs/DMRs annotation
#' 
#' @param tab cpg table
#' @param man manifest
#' 
#' @return plot
#' 
#' @export
#' 
barplots_MM280<-function(tab=NULL,man=NULL){
  
  colors<-RColorBrewer::brewer.pal(8,"Spectral")
  mantmp<-bind_cols(man[,"Feature"],status="MM280",is.sig=FALSE)
  bart<-tab %>% dplyr::select(Feature,status,is.sig) %>% bind_rows(mantmp)
  tmpbart <- bart %>% filter(is.sig) %>% mutate(status=paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) %>% 
    mutate(status=factor(status, levels=c("hypo*", "hypo","plateform","MM280","hyper","hyper*")))
  
  bp1<-bart %>% ggplot(aes(fill=Feature,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  bp2<-bart %>% filter( !grepl("tss",Feature) ) %>%  ggplot(aes(fill=Feature,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  print(bp1)
  print(bp2)
}

#' barplots_450k
#' 
#' Display repartition of DMPs/DMRs annotation
#' 
#' @param tab cpg table
#' @param man manifest
#' 
#' @return plot
#' 
#' @export
#' 
barplots_450k<-function(tab=NULL,man=NULL){
  
  colors<-RColorBrewer::brewer.pal(8,"Spectral")
  mantmp<-bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")],status="Hm450k",is.sig=FALSE)
  bart<-tab %>% dplyr::select(UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island,status,is.sig) %>% bind_rows(mantmp)
  tmpbart <- bart %>% filter(is.sig) %>% mutate(status=paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) %>% 
    mutate(status=factor(status, levels=c("hypo*", "hypo","Hm450k","hyper","hyper*"))) 
  
  bp1<-bart %>% filter(!is.na(status)) %>%
    separate_rows(UCSC_RefGene_Group, sep = ";", convert = FALSE) %>% 
    ggplot(aes(fill=UCSC_RefGene_Group,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  bp2<-bart %>% filter(!is.na(status)) %>%
    separate_rows(Relation_to_UCSC_CpG_Island, sep = ";", convert = FALSE) %>%
    ggplot(aes(fill=Relation_to_UCSC_CpG_Island,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  print(bp1)
  print(bp2)
}


#' barplots_EPIC
#' 
#' Display repartition of DMPs/DMRs annotation
#' 
#' @param tab cpg table
#' @param man manifest
#' 
#' @return plot
#' 
#' @export
#' 
barplots_EPIC<-function(tab=NULL,man=NULL){
  
  cat("\ndebug\n")
  colors<-RColorBrewer::brewer.pal(8,"Spectral")
  mantmp<-bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")],status="EPIC",is.sig=FALSE)
  bart<-tab %>% dplyr::select(UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island,status,is.sig) %>% bind_rows(mantmp)
  tmpbart <- bart %>% filter(is.sig) %>% mutate(status=paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) %>% 
    mutate(status=factor(status, levels=c("hypo*", "hypo","EPIC","hyper","hyper*"))) 
  
  bp1<-bart %>% filter(!is.na(status)) %>%
    separate_rows(UCSC_RefGene_Group, sep = ";", convert = FALSE) %>% 
    ggplot(aes(fill=UCSC_RefGene_Group,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  bp2<-bart %>% filter(!is.na(status)) %>%
    separate_rows(Relation_to_UCSC_CpG_Island, sep = ";", convert = FALSE) %>%
    ggplot(aes(fill=Relation_to_UCSC_CpG_Island,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  print(bp1)
  print(bp2)
}


#' volcano
#' 
#' Volcano plot of DMPs
#' 
#' 
#' @param tab dmps table (dataframe)
#' 
#' @return plot
#' 
#' @export
#'  
volcano<-function(df,  sig=5e-8 ){

  # to make it work with DMRs table.
  df <- df %>% 
    rename_with(~ case_when( . == "HMFDR" ~ "P.Value",  . == "maxdiff" ~ "deltabetas", TRUE ~ .) )
    #rename(any_of(c("P.Value"="HMFDR","deltabetas"="maxdiff")))
  
  
  p <- ggplot( df, aes(x = deltabetas, y = -log(P.Value) ) ) + 
  geom_point( aes(color=ifelse(P.Value>0.05,"lightgray",ifelse(deltabetas>0,"blue","red"))), alpha=0.5 ) +
  scale_colour_manual(values = c("blue", "lightgray", "red")) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed", col="red") +
  guides(color = "none") +
  theme_minimal()

  return(p)
}

#' Manhattan
#' 
#' Manhattan plot of DMPs
#' 
#' The table must contain a column 'chr' for chromosome
#' and a column MAPINFO for the chromosomique position.
#' 
#' @param tab dmps table (dataframe)
#' @param sig cutoff significance
#' 
#' @return plot
#' 
#' @export
#' 
manhattan<-function(df,  sig=5e-8 ){
  
  df <- df  %>% 
    #rename(any_of(c("P.Value"="HMFDR","chr"="seqnames","MAPINFO"="start"))) %>%
    rename_with(~ case_when( . == "HMFDR" ~ "P.Value",  . == "seqnames" ~ "chr", TRUE ~ .) ) %>%
    mutate(MAPINFO = if (exists('MAPINFO', where=.)) MAPINFO else start+(width/2)) %>%
    mutate( chr = factor(chr,levels=c(paste0("chr",seq(1:22)),"chrX","chrY") ) )
  
  data_cum <- df %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(MAPINFO)) %>% 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
    dplyr::select(chr, bp_add)
  
  gwas_data <- df  %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = MAPINFO + bp_add)
  
  axis_set <- gwas_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
  
  ylim <- gwas_data %>% 
    filter(P.Value == min(P.Value)) %>% 
    mutate(ylim = abs(floor(log10(P.Value))) + 2) %>% 
    pull(ylim)
  
  p<-ggplot(gwas_data, aes(x = bp_cum, y = -log10(P.Value), 
                           color = haven::as_factor(chr), size = -log10(P.Value))) +
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(p)") + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      #axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    )
  
  return(p)
}




#######################
#histone marker plot

# mk_hmplot<-function(dm, out, process_id, celllines){
# 
# 	require(data.table)
# 
# ######Extraction infos chrom HMM
# 	dt_chn <- as.data.table(dm)[ grep( "chromatin", annot.type ), .( cpg, annot.type ) ]
# 	dt_chn[ , c( "instance" ) := tstrsplit( annot.type, "_", keep = 3 ) ]
# 	dt_chn[ , c( "cellType", "State" ) := tstrsplit( instance, '-' ) ]
# 
# 	## cf. https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeBroadHmm&hgta_table=wgEncodeBroadHmmGm12878HMM&hgta_doSchema=describe+table+schema
# 
# 	# #### check that one value per cpg per cell line
# 	#if ( dt_chn[ , .N, by = c( 'cpg', 'cellType' ) ][ , max( N ) ] > 1){ stop( "\r Chromatin statuses error!" ) }
# 
# 	v_nb_state <- c( "01_ActivePromoter", "11_Heterochrom/lo", "06_Insulator", "03_PoisedPromoter", "12_Repetitive/CNV", "10_Repressed",
# 		         "04_StrongEnhancer", "08_TxnElongation", "07_TxnTransition", "05_WeakEnhancer", "02_WeakPromoter",  "09_WeakTxn" )
# 	v_names <- c( "Active Promoter", "Weak Promoter", "Inactive/poised Promoter", "Strong enhancer", "Weak/poised enhancer", "Insulator", "Transcriptional transition",
# 		      "Transcriptional elongation", "Weak transcribed", "Polycomb-repressed", "Heterochromatin; low signal", "Repetitive/Copy Number Variation" )
# 	v_cols <- c( "red2", "indianred2", "violetred2", "gold2", "yellow2", "cyan2", "green4", "green4", "palegreen2", "grey35", "grey70", "grey70" )
# 
# 	dt_chn[ , ord_state := factor( State ) ]
# 	levels( dt_chn$ord_state ) <- v_nb_state
# 	dt_chn[ , ord_state := as.factor( as.character( ord_state ) ) ]
# 
# 	dt_sum_chn <- dt_chn[ , .N, by = c( "cellType", "ord_state" ) ]
# 	dt_sum_chn[ , frac := N / sum(N), by = cellType ]
# 	dt_sum_chn<-dt_sum_chn[ cellType %in% celllines ,]
# 	rm( dt_chn )
# 
# 	#chromHMM
# 	p<-ggplot( dt_sum_chn ) +
# 	    geom_col( aes( x = cellType, y = frac * 100, fill = ord_state ), color = "grey40" ) +
# 	    scale_fill_manual( name = "Chromatin Status", labels = v_names, values = v_cols ) +
# 	    xlab( label = "Cell type" ) +
# 	    ylab( label = "Percent" ) +
# 	    theme_minimal() +
# 	    theme( 
# 	    axis.ticks.y = element_line( colour = "black" ), 
# 	    axis.text = element_text( angle=50,hjust=1, size = 12 ),
# 	    axis.title = element_text( size = 15 ) )
# 
# 	return(p)
# 
# }





