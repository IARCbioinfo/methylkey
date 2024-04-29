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
violin_plot<-function(betas, group, numPositions=1000){
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  #cpg<-cbind(group=as.character(group),data.table(t(cpg)) )
  #cpg<-melt( cpg, id.vars="group" )
  #group <- as.character(group)
  cpg <- data.frame(group, t(cpg)) %>% 
    pivot_longer(!group, names_to = "Probe_ID", values_to = "betas")
  
  #colnames(cpg)<-c("group","Probe_ID","betas")
  
  p1<-ggplot(cpg, aes(x=group, y=betas, fill=group)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("")
  
  return(p1)
}

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
barplots_distToSS<-function(dmps, bin=50){
  
  colors<-RColorBrewer::brewer.pal(5,"Spectral")
  
  foo<-dmps %>%
    mutate(status=ifelse(deltabetas > 0, "hyper", "hypo")) %>%
    mutate(status=ifelse(adj.P.Val < 0.05, paste0(status,"*"), status)) %>%
    separate_rows(distToTSS, sep = ";") %>% mutate(distToTSS=as.numeric(distToTSS)) %>%
    mutate(distToTSS=cut_number(distToTSS, bin)) 

  bounds <- stringr::str_match(levels(foo$distToTSS), "\\[?\\(*(.*),(.*)\\]")
  bounds <- data.frame(name=bounds[,1],lower=as.numeric(bounds[,2]), upper=as.numeric(bounds[,3])) %>%
  mutate(labels=case_when(
    lower<0 & upper >0 ~ "TSS",
    lower< -200 & upper > -200 ~ "-200",
    lower< -500 & upper > -500 ~ "-500",
    lower< -1000 & upper > -1000 ~ "-1000",
    lower< 200 & upper > 200 ~ "200",
    lower< 500 & upper > 500 ~ "500",
    lower< 1000 & upper > 1000 ~ "1000",
    TRUE ~ ""
  ))

  moy=table(foo$distToTSS) %>% mean() / 2

  foo %>% mutate(status = factor(status, levels = c("hyper*", "hyper", "hypo*", "hypo")) ) %>%
    group_by(status, distToTSS) %>% count() %>%
    mutate(n=ifelse(grepl(status,"hypo"),-n,n)) %>%
    ggplot(aes(fill=status,x=distToTSS, y=n)) + 
    geom_bar(position="stack", stat="identity") +
    geom_hline(yintercept=moy, color="red") + geom_hline(yintercept=-moy,color="red") +
    labs(x = "distance To Tss", y = "Count", title = "Methylation change by regions") +
    theme_minimal() +
    scale_x_discrete(labels = bounds$labels) +
    theme(legend.text = element_text(hjust = 1,size=12), axis.text.x = element_text(angle=90, hjust=1, size =12) )

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
  mantmp<-bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_Island")],status="Hm450k",is.sig=FALSE)
  bart<-tab %>% dplyr::select(UCSC_RefGene_Group,Relation_to_Island,status,is.sig) %>% bind_rows(mantmp)
  tmpbart <- bart %>% filter(is.sig) %>% mutate(status=paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) %>% 
    mutate(status=factor(status, levels=c("hypo*", "hypo","Hm450k","hyper","hyper*"))) 
  
  bp1<-bart %>% filter(!is.na(status)) %>%
    tidyr::separate_rows(UCSC_RefGene_Group, sep = ";", convert = FALSE) %>% 
    ggplot(aes(fill=UCSC_RefGene_Group,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  bp2<-bart %>% filter(!is.na(status)) %>%
    tidyr::separate_rows(Relation_to_Island, sep = ";", convert = FALSE) %>%
    ggplot(aes(fill=Relation_to_Island,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
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
  mantmp<-bind_cols(man[,c("UCSC_RefGene_Group","Relation_to_Island")],status="EPIC",is.sig=FALSE)
  bart<-tab %>% dplyr::select(UCSC_RefGene_Group,Relation_to_Island,status,is.sig) %>% bind_rows(mantmp)
  tmpbart <- bart %>% filter(is.sig) %>% mutate(status=paste0(status,"*"))
  bart    <- bind_rows(bart, tmpbart ) %>% 
    mutate(status=factor(status, levels=c("hypo*", "hypo","EPIC","hyper","hyper*"))) 
  
  bp1<-bart %>% filter(!is.na(status)) %>%
    separate_rows(UCSC_RefGene_Group, sep = ";", convert = FALSE) %>% 
    ggplot(aes(fill=UCSC_RefGene_Group,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
  bp2<-bart %>% filter(!is.na(status)) %>%
    tidyr::separate_rows(Relation_to_Island, sep = ";", convert = FALSE) %>%
    ggplot(aes(fill=Relation_to_Island,x=status)) + geom_bar(position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme(legend.text=element_text(size=12))
  
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
  df <- df %>% rename_with(~ case_when( . == "HMFDR" ~ "adj.P.Val",  . == "maxdiff" ~ "deltabetas", TRUE ~ .) )

  p <- ggplot( df, aes(x = deltabetas, y = -log10(adj.P.Val) ) ) + 
  geom_point( aes(color=ifelse(adj.P.Val>0.05,"lightgray",ifelse(deltabetas>0,"blue","red"))), alpha=0.5 ) +
  scale_colour_manual(values = c("blue", "lightgray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red") +
  guides(color = "none") +
  theme_minimal()

  return(p)
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
manhattan<-function(df,  sig=5e-8 ){
  
  df <- df  %>%
    #rename(any_of(c("P.Value"="HMFDR","chr"="seqnames","pos"="start"))) %>%
    rename_with(~ case_when( . == "HMFDR" ~ "adj.P.Val",  . == "seqnames" ~ "chr", TRUE ~ .) ) %>%
    mutate(pos = if (exists('pos', where=.)) pos else start+(width/2)) %>%
    mutate( chr = factor(chr,levels=c(paste0("chr",seq(1:22)),"chrX","chrY") ) )
  
  data_cum <- df %>% dplyr::filter(!is.na(chr)) %>%
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>% 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
    dplyr::select(chr, bp_add) 
    
  
  gwas_data <- df  %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = sum(pos,bp_add) )
  
  axis_set <- gwas_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
  
  ylim = min(gwas_data$adj.P.Val)
  sig = abs(floor(log10(ylim))) + 2
  
  p<-ggplot(gwas_data, aes(x = bp_cum, y = -log10(adj.P.Val), 
                           color = haven::as_factor(chr), size = -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, y = "-log<sub>10</sub>(p)") + 
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


#' Create a scatter plot with custom annotations
#'
#' This function generates a scatter plot with custom annotations using ggplot2. The plot displays data from a data frame, where 'mean_oob_grn' and 'mean_oob_red' are used as x and y coordinates, and 'name' is used for annotations.
#'
#' @param df A data frame containing the data to be plotted, with columns 'mean_oob_grn', 'mean_oob_red', and 'name'.
#'
#' @return A ggplot2 scatter plot with custom annotations.
#'
#' @export
plot_background<-function(df){
  ggplot(df,
         aes(x = mean_oob_grn, y= mean_oob_red, label = name)) +
    geom_point() + geom_text(hjust = -0.1, vjust = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
    xlab('Green Background') + ylab('Red Background')
}


#' Create a scatter plot to visualize channel switching
#'
#' This function generates a scatter plot to visualize channel switching using ggplot2. It plots data from a data frame, where 'InfI_switch_G2R' and 'InfI_switch_R2G' are used as x and y coordinates.
#'
#' @param df A data frame containing the data to be plotted, with columns 'InfI_switch_G2R' and 'InfI_switch_R2G'.
#'
#' @return A ggplot2 scatter plot for visualizing channel switching.
#'
#' @export
plot_channel_switch<-function(df){
  ggplot(df) + geom_point(aes(InfI_switch_G2R, InfI_switch_R2G))
}

#' Create a bar plot to visualize missing values (NA)
#'
#' This function generates a bar plot to visualize the number of missing values (NA) for different categories using ggplot2. It plots data from a data frame, where 'name' is used on the x-axis and 'num_na_cg' on the y-axis.
#'
#' @param df A data frame containing the data to be plotted, with columns 'name' and 'num_na_cg'.
#'
#' @return A ggplot2 bar plot for visualizing missing values (NA).
#'
#' @export
plot_NA<-function(df){
  ggplot(df) + geom_bar(aes(x=name,y=num_na_cg), stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' Create a bar plot to visualize detection fractions
#'
#' This function generates a bar plot to visualize detection fractions using ggplot2. It plots data from a data frame, where 'name' is used on the x-axis and 'frac_dt' on the y-axis.
#'
#' @param df A data frame containing the data to be plotted, with columns 'name' and 'frac_dt'.
#'
#' @return A ggplot2 bar plot for visualizing detection fractions.
#'
#' @export
plot_Detection<-function(df){
  ggplot(df) + geom_bar(aes(x=name,y=frac_dt), stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' Create a volcano plot to visualize differential methylation analysis results
#'
#' This function generates a volcano plot to visualize differential methylation analysis results using ggplot2. It takes a data frame 'dt' as input, groups the data by 'dmrdmrtool', 'ID', and 'fdr', calculates summary statistics, and creates a plot.
#'
#' @param df A data frame containing the data for differential methylation analysis.
#'
#' @return A ggplot2 volcano plot showing differential methylation analysis results.
#'
#' @export
my_volcanomean<-function(df){
  volcano <- df %>% group_by(dmrtool,ID,fdr) %>% 
    summarise( no.cpgs=n(), mean.deltabetas=mean(deltabetas), .groups='drop' ) %>%
    ggplot( aes(x = mean.deltabetas, y = -log10(fdr)) ) +
    geom_point( aes(color=ifelse(fdr>0.05,"3notsig",ifelse(mean.deltabetas>0,"1up","2down"))), alpha=0.5 ) +
    scale_colour_manual(values = c("blue", "red", "lightgray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red") +
    guides(color = "none") + facet_grid(. ~ dmrtool) +
    theme_minimal() + ggtitle("A : mean deltabetas")
  return(volcano)
}

#' Create a volcano plot to visualize differential methylation analysis results with harmonized mean
#'
#' This function generates a volcano plot to visualize differential methylation analysis results using ggplot2. It takes a data frame 'dt' as input, groups the data by 'dmrtool', 'ID', and 'fdr', calculates summary statistics with a modified 'mean.deltabetas' calculation, and creates a plot.
#'
#' @param df A data frame containing the data for differential methylation analysis.
#'
#' @return A ggplot2 volcano plot showing differential methylation analysis results with modified mean calculation.
#'
#' @export
my_volcanohmean<-function(df){
  volcano <- df %>% group_by(dmrtool,ID,fdr) %>% 
    summarise( no.cpgs=n(), mean.deltabetas= (n() / sum(1/deltabetas)), .groups='drop' ) %>%
    ggplot( aes(x = mean.deltabetas, y = -log10(fdr)) ) +
    geom_point( aes(color=ifelse(fdr>0.05,"3notsig",ifelse(mean.deltabetas>0,"1up","2down"))), alpha=0.5 ) +
    scale_colour_manual(values = c("blue", "red", "lightgray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red") +
    guides(color = "none") + facet_grid(. ~ dmrtool) +
    theme_minimal() + ggtitle("A : mean deltabetas")
  return(volcano)
}


#' Create a volcano plot to visualize differential methylation analysis results with max value
#'
#' This function generates a volcano plot to visualize differential methylation analysis results using ggplot2. It takes a data frame 'dt' as input, groups the data by 'dmrtool', 'ID', and 'fdr', calculates summary statistics with a modified 'mean.deltabetas' calculation, and creates a plot.
#'
#' @param df A data frame containing the data for differential methylation analysis.
#'
#' @return A ggplot2 volcano plot showing differential methylation analysis results with modified mean calculation.
#'
#' @export
my_volcanomax<-function(df){
  volcano <- df %>% group_by(dmrtool, ID, fdr) %>% 
    summarise( max.deltabetas=max(deltabetas), .groups='drop' ) %>%
    ggplot( aes(x = max.deltabetas, y = -log10(fdr)) ) +
    geom_point( aes(color=ifelse(fdr>0.05,"3notsig",ifelse(max.deltabetas>0,"1up","2down"))), alpha=0.5 ) +
    scale_colour_manual(values = c("blue", "red", "lightgray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red") +
    guides(color = "none") + facet_grid(. ~ dmrtool) +
    theme_minimal() + ggtitle("B : max deltabetas")
  return(volcano)
}

#' Create a volcano plot to visualize statistical results
#'
#' This function generates a volcano plot to visualize statistical results using ggplot2. It takes a data frame 'df' as input, where 'deltabetas' and 'P.Value' are used as x and y coordinates. The points in the plot are colored based on significance and effect size.
#'
#' @param df A data frame containing the data for statistical results with columns 'deltabetas' and 'P.Value'.
#'
#' @return A ggplot2 volcano plot showing statistical results.
#'
#' @export
my_volcano<-function(df){
  
  p <- ggplot( data.frame(df), aes(x = deltabetas, y = -log10(adj.P.Val) ) ) + 
    geom_point( aes(color=ifelse(adj.P.Val>0.05,"lightgray",ifelse(deltabetas>0,"blue","red"))), alpha=0.5 ) +
    scale_colour_manual(values = c("blue", "lightgray", "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col="red") +
    guides(color = "none") +
    theme_minimal()
  
  return(p)
}

#' Create a Venn diagram to visualize intersections of probe IDs
#'
#' This function generates a Venn diagram to visualize the intersections of probe IDs from different sets. It takes a data frame 'dt' as input and a list 'a' containing the names of sets to be intersected.
#'
#' @param dt A data frame containing the data.
#' @param a A list of set names to be intersected.
#'
#' @return A ggvenn Venn diagram showing the intersections between sets.
#'
#' @export
my_venn_<-function(dt,a){
  
  b <- list(`DMRcate` = dt %>% dplyr::filter(dmrtool=="dmrcate") %>% dplyr::pull(Probe_ID) %>% unique(),
            `DMRff` = dt %>% dplyr::filter(dmrtool=="dmrff") %>% dplyr::pull(Probe_ID) %>% unique(),
            `combp` = dt %>% dplyr::filter(dmrtool=="combp") %>% dplyr::pull(Probe_ID) %>% unique(),
            `ipdmr` = dt %>% dplyr::filter(dmrtool=="ipdmr") %>% dplyr::pull(Probe_ID) %>% unique()
  )
  a<-Filter(length, a)
  b<-Filter(length, b)
  venn<-ggvenn::ggvenn( b, names(b) ) + ggtitle("C")
  return(venn)
}

#' Create a Venn diagram to visualize intersections of IDs for different dmrtools
#'
#' This function generates a Venn diagram to visualize the intersections of IDs obtained from different analysis dmrtools. It takes a data frame 'dt' as input and automatically extracts IDs for each dmrtool before passing them to the 'my_venn_' function.
#'
#' @param dt A data frame containing the data.
#'
#' @return A ggvenn Venn diagram showing the intersections between IDs obtained from different dmrtools.
#'
#' @export
my_venn<-function(dt){
  
  a <- list(`DMRcate` = dt %>% dplyr::filter(dmrtool=="dmrcate") %>% dplyr::pull(ID) %>% unique(),
            `DMRff` = dt %>% dplyr::filter(dmrtool=="dmrff") %>% dplyr::pull(ID) %>% unique(),
            `combp` = dt %>% dplyr::filter(dmrtool=="combp") %>% dplyr::pull(ID) %>% unique(),
            `ipdmr` = dt %>% dplyr::filter(dmrtool=="ipdmr") %>% dplyr::pull(ID) %>% unique()
  )
  return(my_venn_(dt,a))
}

#' Create a Venn diagram to visualize intersections of IDs for different dmrtools
#'
#' This function generates a Venn diagram to visualize the intersections of IDs obtained from different analysis dmrtools. It takes a data frame 'dt' as input and automatically extracts IDs for each dmrtool before passing them to the 'my_venn_' function.
#'
#' @param dt A data frame containing the data.
#'
#' @return A ggvenn Venn diagram showing the intersections between IDs obtained from different dmrtools.
#'
#' @export
my_venn_oriented<-function(dt){
  
  a <- list(`DMRcate` = dt %>% dplyr::filter(dmrtool=="dmrcate") %>% dplyr::pull(IDo) %>% unique(),
            `DMRff` = dt %>% dplyr::filter(dmrtool=="dmrff") %>% dplyr::pull(IDo) %>% unique(),
            `combp` = dt %>% dplyr::filter(dmrtool=="combp") %>% dplyr::pull(IDo) %>% unique(),
            `ipdmr` = dt %>% dplyr::filter(dmrtool=="ipdmr") %>% dplyr::pull(IDo) %>% unique()
  )
  return(my_venn_(dt,a))
}

#' Create a density plot of DMR lengths grouped by dmrtool
#'
#' This function generates a density plot of DMR (Differentially Methylated Region) lengths grouped by 'dmrtool'. It takes a data frame 'dt' as input, calculates the length of each DMR, and then creates the density plot.
#'
#' @param dt A data frame containing the data.
#'
#' @return A ggplot2 density plot showing DMR length distributions grouped by 'dmrtool'.
#'
#' @export
my_density<-function(dt){
  
  density <- dt %>% group_by(dmrtool,ID) %>% 
    summarise( no.cpgs=n(), .groups='drop' ) %>%
    separate(ID, into=c("chr","start","end"), sep="-") %>%
    mutate(start=as.numeric(start), end=as.numeric(end)) %>%
    mutate(length=end-start) %>%
    ggdensity(x = "length",
              add = "mean", rug = TRUE,
              color = "dmrtool", fill = "dmrtool",
              palette = c("#00AFBB", "#E7B800", "#FC4E07", "#59981A"),
              title = "D") + scale_x_continuous(trans='log2')
  
  return(density)
}


#' Create a violin plot of DMR lengths grouped by dmrtool
#'
#' This function generates a violin plot of DMR (Differentially Methylated Region) lengths grouped by 'dmrtool'. It takes a data frame 'dt' as input, calculates the length of each DMR, and then creates the violin plot.
#'
#' @param dt A data frame containing the data.
#'
#' @return A ggplot2 violin plot showing DMR length distributions grouped by 'dmrtool'.
#'
#' @export
my_violin<-function(dt){
  violin<-dt %>% group_by(dmrtool,ID) %>% 
    separate(ID, into=c("chr","start","end"), sep="-") %>%
    mutate(start=as.numeric(start), end=as.numeric(end)) %>%
    #summarise( no.cpgs=n(), .groups='drop' ) %>%
    mutate(length=end-start) %>%
    ggviolin(x = "dmrtool", y = "length", fill = "dmrtool",
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#59981A"),
             title = "E")
  return(violin)
}

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
my_track<-function(dt,betas,genome="hg19",chromosome="Chr1",start=0,end=5000){
  
  cpgs <- dt %>% dplyr::pull(Probe_ID) %>% unlist()
  
  foo <- dt %>% group_by(Probe_ID,chr,pos) %>% 
    summarise(.groups='drop') %>% 
    merge(betas[cpgs,] %>% data.frame() %>% 
            rownames_to_column("Probe_ID"), by="Probe_ID") %>% 
    arrange(pos) %>% mutate(start=pos, end=pos) %>% dplyr::select(-c(Probe_ID,pos))
  
  if(nrow(foo)<1) return(NULL)
  
  dmrcate<-dt %>% group_by(dmrtool,chr,start,end) %>% summarise( no.cpgs=n(), .groups='drop' ) %>% filter(dmrtool=="dmrcate")
  dmrff<-dt %>% group_by(dmrtool,chr,start,end) %>% summarise( no.cpgs=n(), .groups='drop' ) %>% filter(dmrtool=="dmrff")
  combp<-dt %>% group_by(dmrtool,chr,start,end) %>% summarise( no.cpgs=n(), .groups='drop' ) %>% filter(dmrtool=="combp")
  ipdmr<-dt %>% group_by(dmrtool,chr,start,end) %>% summarise( no.cpgs=n(), .groups='drop' ) %>% filter(dmrtool=="ipdmr")
  
  #Ideogram track
  itrack <- IdeogramTrack(genome = genome, chromosome = chromosome)
  #Annotation track, title ="CpG"
  atrack0 <- AnnotationTrack(cpgIslands, genome = genome, name = "CpG", chr = chromosome, from = start, to = end)
  #Annotation track
  aTrack1 <- AnnotationTrack(range = dmrcate, name = "DMRcate", chr = chromosome, from = start, to = end)
  aTrack2 <- AnnotationTrack(range = dmrff, name = "DMRff", chr = chromosome, from = start, to = end)
  aTrack3 <- AnnotationTrack(range = combp, name = "combp", chr = chromosome, from = start, to = end)
  aTrack4 <- AnnotationTrack(range = ipdmr, name = "ipdmr", chr = chromosome, from = start, to = end)
  
  foo <- makeGRangesFromDataFrame(foo,keep.extra.columns = TRUE)
  
  dtrack<-DataTrack(foo, name = "probes", groups = sampleSheet$Group, type = "confint", showSampleNames = TRUE, cex.sampleNames = 0.6, genome=genome, chr=chromosome)
  plotTracks(dtrack)
  plotTracks(list(itrack, gtrack,dtrack, atrack0,aTrack1,aTrack2,aTrack3,aTrack4), chr=chromosome ,from = start, to = end)
  
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
my_dmrplot<-function(dt,group,title){
  
  ggplot(dt, aes(x=Probe_ID, y=betas, ymin=(betas-sd/2), ymax=(betas+sd/2) , group=get(group), fill=get(group) ) ) +
    geom_line() + 
    geom_point() +
    geom_ribbon(alpha=0.5) +
    theme(axis.text.x=element_text(angle=50,hjust=1, size=10)) +
    labs(fill=group) +
    ggtitle(title)
  
}

#' Create a Scatter Plot for Visualizing Beta Values
#'
#' This function generates a scatter plot for visualizing beta values using multidimensional scaling (MDS). It takes beta values, group information, sample names, color palette, point size, and other optional parameters for customization.
#'
#' @param betas A matrix of beta values (rows for CpG sites, columns for samples).
#' @param group A vector specifying group labels for samples.
#' @param sample_names A vector of sample names to label the points (optional).
#' @param palette A color palette for group labels.
#' @param size The size of points in the scatter plot.
#' @param showlabels Logical, whether to show sample labels or not.
#' @param ellipse Logical, whether to add ellipses around groups (optional).
#'
#' @return A scatter plot of beta values using multidimensional scaling.
#'
#' @export
my_scatterPlot<-function(betas,group,sample_names,palette,size,showlabels,ellipse){
  
  #get lines in betas without na
  subset <- which(!is.na( rowSums(betas) ))
  subset <- sample(subset,size=100000)
  mds<-t(betas[subset,])  %>% dist() %>% cmdscale() %>% as_tibble()
  colnames(mds) <- c("Dim.1", "Dim.2")
  if(!showlabels) sample_names=NULL
  
  #K-means clustering
  #clust <- kmeans(mds, 3)$cluster %>% as.factor()
  #mds <- mds %>% mutate(cluster = clust)
  mds <- mds %>% mutate(grp=group)
  
  # Plot MDS
  ggscatter(mds, x = "Dim.1", y = "Dim.2", 
            label = sample_names,
            color = "grp",
            palette = palette,
            size = size, 
            ellipse = ellipse,
            #ellipse.type = "convex",
            repel = TRUE,
            max.overlaps=50)
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
manhattan<-function(df,  sig=5e-8 ){
  
  df <- df  %>% data.frame() %>% 
    mutate( chr = factor(chr,levels=c(paste0("chr",seq(1:22)),"chrX","chrY") ) )
  
  data_cum <- df %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>% 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
    dplyr::select(chr, bp_add)
  
  gwas_data <- df  %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(bp_cum = pos + bp_add)
  
  axis_set <- gwas_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(bp_cum))
  
  ylim <- gwas_data %>% 
    filter(P.Value == min(P.Value)) %>% 
    mutate(ylim = abs(floor(log10(P.Value))) + 2) %>% 
    dplyr::pull(ylim)
  
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
  if (genome=="mm10"){ require(BSgenome.Mmusculus.UCSC.mm10) ; species=BSgenome.Mmusculus.UCSC.mm10}
  
  chr.len = seqlengths(species)
  chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
  myIdeo <- GRanges(seqnames = names(chr.len), ranges = IRanges(start = 1, chr.len))
  seqlengths(myIdeo)<-myIdeo@ranges@width
  seqlevels(ranges,pruning.mode="coarse")<-seqlevels(myIdeo)
  seqinfo(ranges)<-seqinfo(myIdeo)
  
  g.po <-ranges[ranges$deltabetas < 0 ]
  g.per<-ranges[ranges$deltabetas > 0 ]
  
  p<- ggbio() + circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
  
  if (length(g.po) > 0){
    values(g.po)$id = "hypo"
    p <- p + circle(g.po, geom = "point", size = 1, aes(x = midpoint, y = "deltabetas", color = id), radius = 19, trackWidth = 20) + scale_colour_manual(values = c("magenta", "green")) 
  }
  if (length(g.per) > 0){	
    values(g.per)$id = "hyper"
    p <- p + circle(g.per, geom = "point", size = 1, aes(x = midpoint, y = "deltabetas", color = id), radius = 41, trackWidth = 20)
  }
  # bug in ggbio waiting for fix
  p<- p + circle(myIdeo, geom = "text", aes(label = seqnames), vjust = 0, radius = 60, trackWidth = 7, size=2)
  return(p)
  
}



#' Generic function for plotting channels
#'
#' This generic function provides a common interface for plotting channels. The
#' behavior of this function depends on the class of the input data.
#'
#' @param x Input data. Depending on its class, behavior varies.
#'
#' @return A plot representing the channels.
#'
#' @export
setGeneric("plot_Channels", function(x) standardGeneric("plot_Channels") )


#' Plot channels from a list
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is a list of channel sets. It extracts and organizes
#' the channels, creates a boxplot, and facets the plot by sentrix.
#'
#' @param x A list where each element represents a set of channels.
#'
#' @return A ggplot2 boxplot representing the channels.
#'
#' @export
setMethod("plot_Channels", signature("list"), definition = function(x) {
  sampNames=names(x)
  if( identical(sampNames,names(x)) ) sampNames<-gsub(".*_R","R",sampNames)
  map( x, ~tibble(UG=.$UG, UR=.$UR) ) %>% bind_rows(.id="name") %>%
    mutate(samp=rep(sampNames,each=nrow(x[[1]]) )) %>%
    tidyr::pivot_longer(-c(name,samp),names_to = "channel", values_to = "Intensity") %>%
    mutate(sentrix=gsub("_.*","",name)) %>%
    ggplot() + geom_boxplot(aes(x=samp,y=Intensity,fill=channel)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) ) +
    scale_fill_manual(values = c("#2ecc71","#c0392b")) +
    scale_y_continuous(trans='log2') +
    ylab("Log2 Intensity") + 
    facet_grid(.~sentrix, scales = "free_x")
})

#' Plot channels from an RGChannelSet object
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates boxplots
#' for red and green channels.
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A ggplot2 plot with separate boxplots for red and green channels.
#'
#' @export
setMethod("plot_Channels", signature("RGChannelSet"), definition = function(x) {
  nsamp=ncol(x)
  ylab<-"log2 intensity of both green and red channel"
  par(xaxt='n')
  boxplot(log2(getRed(x)+1), col = "red", boxwex = 0.25, at= 1:nsamp - 0.175, ylab=ylab, labels=sampleSheet$samples, cex=0.5)
  boxplot(log2(getGreen(x)+1), col = "green", boxwex = 0.25, at= 1:nsamp + 0.175, axis=F , add=T, cex=0.5)
  par(xaxt='s')
  axis(1, at=1:nsamp, labels=SummarizedExperiment::colData(x)$Basename, tick=TRUE, las=2, cex.axis=0.8)
  recordPlot()
})


#' Generic function for plotting channels
#'
#' This generic function provides a common interface for plotting channels. The
#' behavior of this function depends on the class of the input data.
#'
#' @param x Input data. Depending on its class, behavior varies.
#'
#' @return A plot representing the channels.
#'
#' @export
setGeneric("plot_Channels2", function(x) standardGeneric("plot_Channels2") )

#' Call Plot channels from an RGChannelSet object for each sentrix
#'
#' This method of \code{plot_Channels} is specifically for plotting channels
#' when the input data is an RGChannelSet object. It creates separate boxplots
#' for each sentrix
#'
#' @param x An RGChannelSet object containing microarray channel data.
#'
#' @return A list of boxplots.
#'
#' @export
setMethod("plot_Channels2", signature("RGChannelSet"), definition = function(x) {
  
  nsubplots <- length(sentrix)
  stx<-gsub("_.*","",colnames(x)) %>% unique()
  
  lapply(1:nsubplots, function(i){
    subsel <- grepl(stx[i], colnames(x))
    substx <- x[, subsel]
    plot_Channels(substx)
  })
  
})


