#require(tidyverse)
#require(ggsankeyfier)
#require(patchwork) 


#' Create a Sankey Plot with Colored Nodes
#'
#' This function generates a Sankey plot from a given dataframe, visualizing relationships 
#' between genes and terms. The node colors are determined by the provided palette.
#'
#' @param df A data frame containing gene-term associations with the following required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'   - `status`: Character, indicating if a node is "up", "down", or another category.
#' @param palette A named vector of colors. The colors should be named according to the 
#'   categories in the `status` column, with special colors for "up" and "down".
#'
#' @return A ggplot object representing the Sankey diagram.
#'
sankeyPlot_<-function(df,palette){
  
  df <- df %>% 
    mutate(logP = -log10(P.value),Term = fct_reorder(Term, gene_count)) %>%
    arrange(gene_count,Term) %>%
    mutate(yy=cumsum(gene_count)) %>% 
    mutate(term_rank = as.numeric(Term)-1, offset = term_rank * max(gene_count) * 2) %>%
    group_by(Term) %>% 
    mutate(yy=mean(yy)+offset ) %>%
    pivot_stages_longer(stages_from = c("Genes", "Term"), values_from = "gene_count",additional_aes_from =c("Term","Odds.Ratio","logP","yy","status")) %>%
    mutate(status=ifelse(stage=="Genes",status,as.character(node)))

  pos <- position_sankey( nudge_x = -0.06, v_space = "auto", order = "ascending")

  sankey_plot <- ggplot(df, aes(x = stage, y = gene_count, group = node, connector = connector, edge_id = edge_id)) +
    geom_sankeynode(v_space = "auto" , aes(fill = status)) +
    geom_sankeyedge(v_space = "auto", aes(fill = Term, label=status)) +
    geom_text( aes(label = node, cex=stage), position = pos, stat = "sankeynode", hjust=1 ) + 
    theme_minimal() +
    scale_x_discrete(expand = expansion(add = c(0.2, 0))) +
    scale_size_manual(values = c(1.5,3),guide="none") +
    scale_fill_manual(values=palette,guide="none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(sankey_plot)
}

#' Create a Dot Plot for Gene-Term Associations
#'
#' This function generates a dot plot that visualizes gene-term associations 
#' based on odds ratio, gene count, and statistical significance.
#'
#' @param df A data frame containing gene-term associations with the following required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'
#' @return A ggplot object representing the dot plot.
#' 
dotPlot_<-function(df){
  
  dotplot_data <- df %>%
    mutate(logP = -log10(P.value),Term = fct_reorder(Term, gene_count)) %>%
    arrange(gene_count,Term) %>%
    mutate(yy=cumsum(gene_count)) %>% 
    mutate(term_rank = as.numeric(Term)-1, offset = term_rank * max(gene_count) * 2) %>%
    group_by(Term) %>% 
    mutate(yy=mean(yy)+offset )
  
  dotplot <- ggplot(dotplot_data, aes(x =Odds.Ratio, y = yy, size = gene_count, color = logP)) +
    geom_point(alpha = 0.8) +
    geom_rect(aes(xmin = min(Odds.Ratio)-0.05*max(Odds.Ratio), ymin = 0, xmax=max(Odds.Ratio)+0.05*max(Odds.Ratio), ymax=max(yy)+0.2*max(yy)), color = "black", fill = NA, linewidth = 0.1 ) +
    scale_color_viridis_c() +
    coord_cartesian(clip="off") +
    theme_minimal() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 9),  
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.y = unit(0.2, "cm"),
          legend.box = "horizontal",
          legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
          legend.position=c(.5,.8)
    )
  
  return(dotplot)
}
  
#' Create a Combined Sankey and Dot Plot
#'
#' This function generates a composite visualization that includes both a 
#' Sankey plot and a dot plot, allowing for an intuitive representation of 
#' gene-term associations.
#'
#' @param df A data frame containing gene-term associations with the following required columns:
#'   - `P.value`: Numeric, p-values for the terms.
#'   - `gene_count`: Numeric, count of genes associated with each term.
#'   - `Term`: Character, names of the terms.
#'   - `Odds.Ratio`: Numeric, odds ratio values.
#'   - `status`: Character, indicating node status (e.g., "up", "down", or other terms).
#' @param palette A named vector of colors for the terms and statuses. If `NULL`, 
#'   the function generates a default palette using `viridis::viridis()`, 
#'   with fixed colors for `"up"` (`#317AC1`) and `"down"` (`#E1A624`).
#'
#' @return A ggplot object combining a Sankey plot and a dot plot.
#' 
#' @export
#' 
sankeyDotPlot<-function(df,palette=NULL){

  if(is.null(palette)){
    palette <- c(
      set_names(viridis::viridis(n_distinct(df$Term)), 
                setdiff(unique(df$Term), c("up", "down"))),
      "up" = "#317AC1",
      "down" = "#E1A624"
    )
  }
  
  sp <- sankeyPlot_(df,palette)
  dp <- dotPlot_(df)
  sdp<-sp + dp + plot_layout(widths = c(5, 3)) & scale_y_continuous(limits = c(0, 3100)) 
  return(sdp)
}


