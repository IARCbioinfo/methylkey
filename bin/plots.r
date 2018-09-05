#####################################################################################
#
# Title  : plots.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 12/07/2018
#
#####################################################################################

#######################
# multiplot
multiplot<-function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#######################
# violin plot
violin_plot<-function(pdata,betas,group,samples,platform, path, out="./",genome="hg19",regions){

	suppressPackageStartupMessages( require(ggplot2) )
	suppressPackageStartupMessages( require(annotatr) )
	require(data.table)

	if (!is.data.table(pdata)){ pdata<-data.table(pdata) }
	#colnames(pdata)[barcode]<-"barcode"
	pdata$barcode<-as.character(pdata$barcode)

	#get cpgs annotation
	if (platform=="IlluminaHumanMethylation450k"){regions=paste0(datadir,"/illumina450k.bed")}
	if (platform=="IlluminaHumanMethylationEPIC"){regions=paste0(datadir,"/illuminaEpic.bed")}
	cpg_annotated<-mk_annotate(regions,genome=genome)

	#format annotation
	annot<-data.table(cbind(cpg_annotated$name, cpg_annotated$annot$symbol, cpg_annotated$annot$type))
	annot<-annot[!grepl("chromatin", V3) ]
	annot<-unique(annot)
	annot[, c("ref","type","annot") := tstrsplit(V3,"_") ]

	#merge betas and pdata in a data.table
	colnames(betas)<-gsub("X","",colnames(betas))
	b<-data.table(melt(as.matrix(betas)))
	b<-merge(pdata[,.( get(samples), grp=get(group))],b,by.x="V1",by.y="Var2")
	if (!is.character(b$Var1)) { b$Var1<-as.character(b$Var1) }
	
	#get cpg islands annotation
	cpg<-merge(b, annot[type=="cpg"], by.x="Var1", by.y="V1")
	#get genes structure annotation
	genes<-merge(b, annot[type=="genes"], by.x="Var1", by.y="V1",allow.cartesian=TRUE)
	#violin plots
	p1<-ggplot(cpg, aes_string(x="grp", y="value", fill="grp")) + geom_violin() + theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + guides(fill =FALSE) +
        ylab("Methylation %") + xlab("")
	p2<-ggplot(cpg, aes_string(x="annot", y="value", fill="grp")) + geom_violin() + theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + guides(fill =FALSE) +
        ylab("Methylation %") + xlab("")
	p3<-ggplot(genes, aes_string(x="annot", y="value", fill="grp")) + geom_violin(scale="area") + theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) +
        ylab("Methylation %") + xlab("")

	jpeg(paste0(out, "/violin1.jpg"), width=1000, height=800)
	print(p1)
	dev.off()

	jpeg(paste0(out, "/violin2.jpg"), width=1000, height=800)
	print(p2)
	dev.off()

	jpeg(paste0(out, "/violin3.jpg"), width=1000, height=800)
	print(p3)
	dev.off()

	#jpeg(paste(out, "violin.jpg", sep="/"), width=1000, height=800)
	#multiplot(p1,p2,p3,cols=3)
	#dev.off()
}

#######################
#plot a circus manhattan plot
circusplot<-function(ranges, genome){

	suppressPackageStartupMessages( require(ggbio) )

	if (genome=="hg38"){ require(BSgenome.Hsapiens.UCSC.hg38) ; species=BSgenome.Hsapiens.UCSC.hg38}
	if (genome=="hg19"){ require(BSgenome.Hsapiens.UCSC.hg19) ; species=BSgenome.Hsapiens.UCSC.hg19}

	chr.len = seqlengths(species)
	chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
	myIdeo <- GRanges(seqnames = names(chr.len), ranges = IRanges(start = 1, chr.len))
	seqlengths(myIdeo)<-myIdeo@ranges@width

	seqlevels(ranges) = names(chr.len)
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





#######################
#mk_barplot

mk_barplot<-function(annotated_regions, betafc, what=c("cpgi","genes") ){
	
	suppressPackageStartupMessages( require(ggplot2) )
	require(RColorBrewer)

	annotated_regions$DM_Status[ betafc > 0 ] = "hyper"
	annotated_regions$DM_Status[ betafc < 0 ] = "hypo"

	if (what=="cpgi"){
		annots_order = paste0( genome, c('_cpg_southShelves', '_cpg_southShores', '_cpg_islands', '_cpg_northShores', '_cpg_northShelves', '_cpg_inter'))
		x_order = c('hyper','hypo')
		colors<-brewer.pal(12,"Paired")[c(1,7,8,7,1,2)]
	}
	else {
		annots_order = paste0( genome, c('_genes_intergenic', '_genes_1to5kb', '_genes_5UTRs', '_genes_promoters', '_genes_firstexons', '_genes_exons', '_genes_introns'))
		x_order = c('hyper','hypo')
		colors<-brewer.pal(7, "Paired")
	}

	foo<-annotated_regions[ annotated_regions$annot.type %in% annots_order, ]
	foo<-rbind(data.frame(type=as.factor(foo$annot.type), status=foo$DM_Status, val=1), data.frame(type=foo$annot.type, status="all", val=1) )
	levels(foo$type)<-annots_order
	
	p<-ggplot(foo, aes(fill=type, x=status, y=val)) + geom_bar( stat="identity", position="fill") + coord_flip() + scale_fill_manual(values=colors) + theme		(legend.text=element_text(size=12))
	return(p)
}


#######################
#histone marker plot

mk_hmplot<-function(dm, out, process_id, celllines){

	require(data.table)

######Extraction infos chrom HMM
	dt_chn <- as.data.table(dm)[ grep( "chromatin", annot.type ), .( cpg, annot.type ) ]
	dt_chn[ , c( "instance" ) := tstrsplit( annot.type, "_", keep = 3 ) ]
	dt_chn[ , c( "cellType", "State" ) := tstrsplit( instance, '-' ) ]

	## cf. https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeBroadHmm&hgta_table=wgEncodeBroadHmmGm12878HMM&hgta_doSchema=describe+table+schema

	# #### check that one value per cpg per cell line
	#if ( dt_chn[ , .N, by = c( 'cpg', 'cellType' ) ][ , max( N ) ] > 1){ stop( "\r Chromatin statuses error!" ) }

	v_nb_state <- c( "01_ActivePromoter", "11_Heterochrom/lo", "06_Insulator", "03_PoisedPromoter", "12_Repetitive/CNV", "10_Repressed",
		         "04_StrongEnhancer", "08_TxnElongation", "07_TxnTransition", "05_WeakEnhancer", "02_WeakPromoter",  "09_WeakTxn" )
	v_names <- c( "Active Promoter", "Weak Promoter", "Inactive/poised Promoter", "Strong enhancer", "Weak/poised enhancer", "Insulator", "Transcriptional transition",
		      "Transcriptional elongation", "Weak transcribed", "Polycomb-repressed", "Heterochromatin; low signal", "Repetitive/Copy Number Variation" )
	v_cols <- c( "red2", "indianred2", "violetred2", "gold2", "yellow2", "cyan2", "green4", "green4", "palegreen2", "grey35", "grey70", "grey70" )

	dt_chn[ , ord_state := factor( State ) ]
	levels( dt_chn$ord_state ) <- v_nb_state
	dt_chn[ , ord_state := as.factor( as.character( ord_state ) ) ]

	dt_sum_chn <- dt_chn[ , .N, by = c( "cellType", "ord_state" ) ]
	dt_sum_chn[ , frac := N / sum(N), by = cellType ]
	dt_sum_chn<-dt_sum_chn[ cellType %in% celllines ,]
	rm( dt_chn )

	#chromHMM
	p<-ggplot( dt_sum_chn ) +
	    geom_col( aes( x = cellType, y = frac * 100, fill = ord_state ), color = "grey40" ) +
	    scale_fill_manual( name = "Chromatin Status", labels = v_names, values = v_cols ) +
	    xlab( label = "Cell type" ) +
	    ylab( label = "Percent" ) +
	    theme_minimal() +
	    theme( 
	    axis.ticks.y = element_line( colour = "black" ), 
	    axis.text = element_text( angle=50,hjust=1, size = 12 ),
	    axis.title = element_text( size = 15 ) )

	return(p)

}



