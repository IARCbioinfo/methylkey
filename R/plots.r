#' barplots_distToSS
#' 
#' Display barplot of deltabetas distribution around TSS
#' 
#' @param dmps data.frame of dmps
#' @param bin size of bin
#' 
#' @return plot
#' 
#' @export
#' 
barplots_distToSS <- function(dmps, bin = 50){
  
  colors <- RColorBrewer::brewer.pal(5,"Spectral")
  
  foo <- dmps |>
    mutate(status = ifelse(deltabetas > 0, "hyper", "hypo")) |>
    mutate(status = ifelse(adj.P.Val < 0.05, paste0(status,"*"), status)) |>
    separate_rows(distToTSS, sep  = ";") |> mutate(distToTSS = as.numeric(distToTSS)) |>
    mutate(distToTSS = cut_number(distToTSS, bin)) 

  bounds <- stringr::str_match(levels(foo$distToTSS), "\\[?\\(*(.*),(.*)\\]")
  bounds <- data.frame(name = bounds[,1],lower = as.numeric(bounds[,2]), upper = as.numeric(bounds[,3])) |>
  mutate(labels = case_when(
    lower<0 & upper >0 ~ "TSS",
    lower< -200 & upper > -200 ~ "-200",
    lower< -500 & upper > -500 ~ "-500",
    lower< -1000 & upper > -1000 ~ "-1000",
    lower< 200 & upper > 200 ~ "200",
    lower< 500 & upper > 500 ~ "500",
    lower< 1000 & upper > 1000 ~ "1000",
    TRUE ~ ""
  ))

  moy = table(foo$distToTSS) |> mean() / 2

  foo |> mutate(status  = factor(status, levels  = c("hyper*", "hyper", "hypo*", "hypo")) ) |>
    group_by(status, distToTSS) |> count() |>
    mutate(n = ifelse(grepl(status,"hypo"),-n,n)) |>
    ggplot(aes(fill = status,x = distToTSS, y = n)) + 
    geom_bar(position = "stack", stat = "identity") +
    geom_hline(yintercept = moy, color = "red") + geom_hline(yintercept = -moy,color = "red") +
    labs(x  = "distance To Tss", y  = "Count", title  = "Methylation change by regions") +
    theme_minimal() +
    scale_x_discrete(labels  = bounds$labels) +
    theme(legend.text  = element_text(hjust  = 1,size = 12), axis.text.x  = element_text(angle = 90, hjust = 1, size  = 12) )

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
barplots_ <- function(plateform = "IlluminaHumanMethylation450k",tab = NULL,man = NULL){
  
  cat(paste0("\nLocation annotation for ",plateform," probes, hyper/hypo methylated probes and significantly hyper/hypo methylated probes\n"))
  if( plateform== "IlluminaHumanMethylation450k") { barplots_450k(tab,man) }
  if( plateform== "IlluminaMouseMethylation285k") { barplots_MM280(tab,man) }
  if( plateform== "IlluminaHumanMethylationEPIC") { barplots_EPIC(tab,man) }
                   
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
barplots_MM280 <- function(tab = NULL,man = NULL){
  
  colors <- RColorBrewer::brewer.pal(8,"Spectral")
  mantmp <- bind_cols(man[,"Feature"],status = "MM280",is.sig = FALSE)
  bart <- tab |> dplyr::select(Feature,status,is.sig) |> bind_rows(mantmp)
  tmpbart <- bart |> filter(is.sig) |> mutate(status = paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) |> 
    mutate(status = factor(status, levels = c("hypo*", "hypo","plateform","MM280","hyper","hyper*")))
  
  bp1 <- bart |> ggplot(aes(fill = Feature,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  bp2 <- bart |> filter( !grepl("tss",Feature) ) |>  ggplot(aes(fill = Feature,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  
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
barplots_450k <- function(tab = NULL,man = NULL){
  
  colors <- RColorBrewer::brewer.pal(8,"Spectral")
  mantmp <- bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_Island")],status = "Hm450k",is.sig = FALSE)
  bart <- tab |> dplyr::select(UCSC_RefGene_Group,Relation_to_Island,status,is.sig) |> bind_rows(mantmp)
  tmpbart <- bart |> filter(is.sig) |> mutate(status = paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) |> 
    mutate(status = factor(status, levels = c("hypo*", "hypo","Hm450k","hyper","hyper*"))) 
  
  bp1 <- bart |> filter(!is.na(status)) |>
    tidyr::separate_rows(UCSC_RefGene_Group, sep  = ";", convert  = FALSE) |> 
    ggplot(aes(fill = UCSC_RefGene_Group,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  
  bp2 <- bart |> filter(!is.na(status)) |>
    tidyr::separate_rows(Relation_to_Island, sep  = ";", convert  = FALSE) |>
    ggplot(aes(fill = Relation_to_Island,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  
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
barplots_EPIC <- function(tab = NULL,man = NULL){
  
  cat("\ndebug\n")
  colors <- RColorBrewer::brewer.pal(8,"Spectral")
  mantmp <- bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_Island")],status = "EPIC",is.sig = FALSE)
  bart <- tab |> dplyr::select(UCSC_RefGene_Group,Relation_to_Island,status,is.sig) |> bind_rows(mantmp)
  tmpbart <- bart |> filter(is.sig) |> mutate(status = paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) |> 
    mutate(status = factor(status, levels = c("hypo*", "hypo","EPIC","hyper","hyper*"))) 
  
  bp1 <- bart |> filter(!is.na(status)) |>
    separate_rows(UCSC_RefGene_Group, sep  = ";", convert  = FALSE) |> 
    ggplot(aes(fill = UCSC_RefGene_Group,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  
  bp2 <- bart |> filter(!is.na(status)) |>
    tidyr::separate_rows(Relation_to_Island, sep  = ";", convert  = FALSE) |>
    ggplot(aes(fill = Relation_to_Island,x = status)) + geom_bar(position = "fill") + coord_flip() + scale_fill_manual(values = colors) + theme(legend.text = element_text(size = 12))
  
  print(bp1)
  print(bp2)
}




#' Manhattan
#' 
#' Manhattan plot of DMPs
#' 
#' The table must contain a column 'chr' for chromosome
#' and a column pos for the chromosomique position.
#' 
#' @param tab dmps table (dataframe)
#' @param sig cutoff significance
#' 
#' @return plot
#' 
#' @export
#' 
manhattan <- function(df,  sig = 5e-8 ){
  
  df <- df  |>
    #rename(any_of(c("P.Value" = "HMFDR","chr" = "seqnames","pos" = "start"))) |>
    rename_with(~ case_when( .== "HMFDR" ~ "adj.P.Val",  .== "seqnames" ~ "chr", TRUE ~ .) ) |>
    mutate(pos  = if (exists('pos', where = .)) pos else start+(width/2)) |>
    mutate( chr  = factor(chr,levels = c(paste0("chr",seq(1:22)),"chrX","chrY") ) )
  
  data_cum <- df |> dplyr::filter(!is.na(chr)) |>
    group_by(chr) |> 
    summarise(max_bp  = max(pos)) |> 
    mutate(bp_add  = lag(cumsum(as.numeric(max_bp)), default  = 0)) |> 
    dplyr::select(chr, bp_add) 
    
  
  gwas_data <- df  |> 
    inner_join(data_cum, by  = "chr") |> 
    mutate(bp_cum  = sum(pos,bp_add) )
  
  axis_set <- gwas_data |> 
    group_by(chr) |> 
    summarize(center  = mean(bp_cum))
  
  ylim  = min(gwas_data$adj.P.Val)
  sig  = abs(floor(log10(ylim))) + 2
  
  p <- ggplot(gwas_data, aes(x  = bp_cum, y  = -log10(adj.P.Val), 
                           color  = haven::as_factor(chr), size  = -log10(adj.P.Val))) +
    geom_hline(yintercept  = -log10(sig), color  = "grey40", linetype  = "dashed") + 
    geom_point(alpha  = 0.75) +
    scale_x_continuous(label  = axis_set$chr, breaks  = axis_set$center) +
    scale_y_continuous(expand  = c(0,0), limits  = c(0, ylim)) +
    scale_color_manual(values  = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
    scale_size_continuous(range  = c(0.5,3)) +
    labs(x  = NULL, y  = "-log<sub>10</sub>(p)") + 
    theme_minimal() +
    theme( 
      legend.position  = "none",
      panel.border  = element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor.x  = element_blank(),
      #axis.title.y  = element_markdown(),
      axis.text.x  = element_text(angle  = 60, size  = 8, vjust  = 0.5)
    )
  
  return(p)
}




#######################
#histone marker plot

# mk_hmplot <- function(dm, out, process_id, celllines){
# 
# 	require(data.table)
# 
# ######Extraction infos chrom HMM
# 	dt_chn <- as.data.table(dm)[ grep( "chromatin", annot.type ), .( cpg, annot.type ) ]
# 	dt_chn[ , c( "instance" ) : = tstrsplit( annot.type, "_", keep  = 3 ) ]
# 	dt_chn[ , c( "cellType", "State" ) : = tstrsplit( instance, '-' ) ]
# 
# 	## cf. https://genome.ucsc.edu/cgi-bin/hgTables?db = hg19&hgta_group = regulation&hgta_track = wgEncodeBroadHmm&hgta_table = wgEncodeBroadHmmGm12878HMM&hgta_doSchema = describe+table+schema
# 
# 	# #### check that one value per cpg per cell line
# 	#if ( dt_chn[ , .N, by  = c( 'cpg', 'cellType' ) ][ , max( N ) ] > 1){ stop( "\r Chromatin statuses error!" ) }
# 
# 	v_nb_state <- c( "01_ActivePromoter", "11_Heterochrom/lo", "06_Insulator", "03_PoisedPromoter", "12_Repetitive/CNV", "10_Repressed",
# 		         "04_StrongEnhancer", "08_TxnElongation", "07_TxnTransition", "05_WeakEnhancer", "02_WeakPromoter",  "09_WeakTxn" )
# 	v_names <- c( "Active Promoter", "Weak Promoter", "Inactive/poised Promoter", "Strong enhancer", "Weak/poised enhancer", "Insulator", "Transcriptional transition",
# 		      "Transcriptional elongation", "Weak transcribed", "Polycomb-repressed", "Heterochromatin; low signal", "Repetitive/Copy Number Variation" )
# 	v_cols <- c( "red2", "indianred2", "violetred2", "gold2", "yellow2", "cyan2", "green4", "green4", "palegreen2", "grey35", "grey70", "grey70" )
# 
# 	dt_chn[ , ord_state : = factor( State ) ]
# 	levels( dt_chn$ord_state ) <- v_nb_state
# 	dt_chn[ , ord_state : = as.factor( as.character( ord_state ) ) ]
# 
# 	dt_sum_chn <- dt_chn[ , .N, by  = c( "cellType", "ord_state" ) ]
# 	dt_sum_chn[ , frac : = N / sum(N), by  = cellType ]
# 	dt_sum_chn <- dt_sum_chn[ cellType %in% celllines ,]
# 	rm( dt_chn )
# 
# 	#chromHMM
# 	p <- ggplot( dt_sum_chn ) +
# 	    geom_col( aes( x  = cellType, y  = frac * 100, fill  = ord_state ), color  = "grey40" ) +
# 	    scale_fill_manual( name  = "Chromatin Status", labels  = v_names, values  = v_cols ) +
# 	    xlab( label  = "Cell type" ) +
# 	    ylab( label  = "Percent" ) +
# 	    theme_minimal() +
# 	    theme( 
# 	    axis.ticks.y  = element_line( colour  = "black" ), 
# 	    axis.text  = element_text( angle = 50,hjust = 1, size  = 12 ),
# 	    axis.title  = element_text( size  = 15 ) )
# 
# 	return(p)
# 
# }












#' Create a genome track plot with CpG, DMRs, and probe annotations
#'
#' This function generates a genome track plot that includes various annotations, such as CpG islands, DMRs (Differentially Methylated Regions), and probe annotations. It takes a data frame 'dt' containing the probe information, a matrix 'betas' containing beta values, and additional parameters for specifying the genome region to display.
#'
#' @param dt A data frame containing the probe information.
#' @param betas A matrix containing beta values.
#' @param genome The genome assembly version (e.g., "hg19").
#' @param chromosome The chromosome to display (e.g., "Chr1").
#' @param start The start position of the genome region to display.
#' @param end The end position of the genome region to display.
#'
#' @return A genome track plot displaying CpG islands, DMRs, and probe annotations.
#'
#' @export
my_track <- function(dt,betas,genome = "hg19",chromosome = "Chr1",start = 0,end = 5000){
  
  cpgs <- dt |> dplyr::pull(Probe_ID) |> unlist()
  
  foo <- dt |> group_by(Probe_ID,chr,pos) |> 
    summarise(.groups = 'drop') |> 
    merge(betas[cpgs,] |> data.frame() |> 
            rownames_to_column("Probe_ID"), by = "Probe_ID") |> 
    arrange(pos) |> mutate(start = pos, end = pos) |> dplyr::select(-c(Probe_ID,pos))
  
  if(nrow(foo)<1) return(NULL)
  
  dmrcate <- dt |> group_by(dmrtool,chr,start,end) |> summarise( no.cpgs = n(), .groups = 'drop' ) |> filter(dmrtool== "dmrcate")
  dmrff <- dt |> group_by(dmrtool,chr,start,end) |> summarise( no.cpgs = n(), .groups = 'drop' ) |> filter(dmrtool== "dmrff")
  combp <- dt |> group_by(dmrtool,chr,start,end) |> summarise( no.cpgs = n(), .groups = 'drop' ) |> filter(dmrtool== "combp")
  ipdmr <- dt |> group_by(dmrtool,chr,start,end) |> summarise( no.cpgs = n(), .groups = 'drop' ) |> filter(dmrtool== "ipdmr")
  
  #Ideogram track
  itrack <- IdeogramTrack(genome  = genome, chromosome  = chromosome)
  #Annotation track, title  = "CpG"
  atrack0 <- AnnotationTrack(cpgIslands, genome  = genome, name  = "CpG", chr  = chromosome, from  = start, to  = end)
  #Annotation track
  aTrack1 <- AnnotationTrack(range  = dmrcate, name  = "DMRcate", chr  = chromosome, from  = start, to  = end)
  aTrack2 <- AnnotationTrack(range  = dmrff, name  = "DMRff", chr  = chromosome, from  = start, to  = end)
  aTrack3 <- AnnotationTrack(range  = combp, name  = "combp", chr  = chromosome, from  = start, to  = end)
  aTrack4 <- AnnotationTrack(range  = ipdmr, name  = "ipdmr", chr  = chromosome, from  = start, to  = end)
  
  foo <- makeGRangesFromDataFrame(foo,keep.extra.columns  = TRUE)
  
  dtrack <- DataTrack(foo, name  = "probes", groups  = sampleSheet$Group, type  = "confint", showSampleNames  = TRUE, cex.sampleNames  = 0.6, genome = genome, chr = chromosome)
  plotTracks(dtrack)
  plotTracks(list(itrack, gtrack,dtrack, atrack0,aTrack1,aTrack2,aTrack3,aTrack4), chr = chromosome ,from  = start, to  = end)
  
}

#' Create a DMR plot for visualizing beta values with confidence intervals
#'
#' This function generates a DMR (Differentially Methylated Region) plot for visualizing beta values with confidence intervals. It takes a data frame 'dt' containing the required data, a grouping variable 'group', and a title for the plot.
#'
#' @param dt A data frame containing the required data.
#' @param group A grouping variable for distinguishing different groups in the plot.
#' @param title A title for the plot.
#'
#' @return A ggplot2 plot displaying beta values with confidence intervals for different groups.
#'
#' @export
my_dmrplot <- function(dt,group,title){
  
  ggplot(dt, aes(x = Probe_ID, y = betas, ymin = (betas-sd/2), ymax = (betas+sd/2) , group = get(group), fill = get(group) ) ) +
    geom_line() + 
    geom_point() +
    geom_ribbon(alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 50,hjust = 1, size = 10)) +
    labs(fill = group) +
    ggtitle(title)
  
}




#' Create a Manhattan Plot for GWAS Data
#'
#' This function generates a Manhattan plot for visualizing genome-wide association study (GWAS) data. It takes a data frame containing chromosome (chr), position (pos), and p-value (P.Value) columns as input.
#'
#' @param df A data frame containing GWAS data with columns chr (chromosome), pos (position), and P.Value (p-value).
#' @param sig The significance threshold for highlighting points on the plot (default is 5e-8).
#'
#' @return A Manhattan plot visualizing GWAS data.
#'
#' @export
manhattan <- function(df,  sig = 5e-8 ){
  
  df <- df  |> data.frame() |> 
    mutate( chr  = factor(chr,levels = c(paste0("chr",seq(1:22)),"chrX","chrY") ) )
  
  data_cum <- df |> 
    group_by(chr) |> 
    summarise(max_bp  = max(pos)) |> 
    mutate(bp_add  = lag(cumsum(as.numeric(max_bp)), default  = 0)) |> 
    dplyr::select(chr, bp_add)
  
  gwas_data <- df  |> 
    inner_join(data_cum, by  = "chr") |> 
    mutate(bp_cum  = pos + bp_add)
  
  axis_set <- gwas_data |> 
    group_by(chr) |> 
    summarize(center  = mean(bp_cum))
  
  ylim <- gwas_data |> 
    filter(P.Value== min(P.Value)) |> 
    mutate(ylim  = abs(floor(log10(P.Value))) + 2) |> 
    dplyr::pull(ylim)
  
  p <- ggplot(gwas_data, aes(x  = bp_cum, y  = -log10(P.Value), 
                           color  = haven::as_factor(chr), size  = -log10(P.Value))) +
    geom_hline(yintercept  = -log10(sig), color  = "grey40", linetype  = "dashed") + 
    geom_point(alpha  = 0.75) +
    scale_x_continuous(label  = axis_set$chr, breaks  = axis_set$center) +
    scale_y_continuous(expand  = c(0,0), limits  = c(0, ylim)) +
    scale_color_manual(values  = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
    scale_size_continuous(range  = c(0.5,3)) +
    labs(x  = NULL, 
         y  = "-log<sub>10</sub>(p)") + 
    theme_minimal() +
    theme( 
      legend.position  = "none",
      panel.border  = element_blank(),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor.x  = element_blank(),
      #axis.title.y  = element_markdown(),
      axis.text.x  = element_text(angle  = 60, size  = 8, vjust  = 0.5)
    )
  
  return(p)
}







