#' PCA analysis
#' 
#' Estimate correlation covariates vs. PCA axes
#' 
#' @param betas matrix of betas
#' @param pdata sampleSheet of data
#' @param nPC maximum number of principal component to display
#' 
#' @return betas matrix
#' 
#' @export
#' 
makepca <- function( betas, pdata, nPC=ncol(betas) ){
  
  nPC=min( nPC,ncol(betas) )
  
  invtbetas = 1/t(betas) # n x p required for prcomp
  #mval<-t(beta2m(as.matrix(betas)))
  invtbetas[!is.finite(invtbetas)]<-min(invtbetas[is.finite(invtbetas)])
  
  #variance for a probe should be > 0
  sel<-which(apply(invtbetas, 2, var)==0)
  
  if (length(sel)) { 
    pca = prcomp(invtbetas[,-sel],scale=T, center=T ) 
  } else { 
    pca = prcomp(invtbetas,scale=T, center=T) 
  }
  
  # plot proportion of var explained
  propVarExpl = summary(pca)$importance[2,]
  
  ###########################
  #plot contributions
  plot_PCA_contribution <- function( pca, nPC = 4, title = "" ){
    v_var <- summary( pca )$importance[ 2, 1:nPC ] 
    return( ggplot2::ggplot(data.frame( PCA_axis = names( v_var ), val = v_var )) 
            + ggplot2::geom_col(aes( x=reorder( PCA_axis, 1:nPC ), y=val )) 
            + labs( title = title, subtitle = "PCA contribution", x = "PC", y = "Prption expl. var." )
            + theme( plot.title = element_text(lineheight=.8, face="bold") ) 
    )
  }
  
  p<-plot_PCA_contribution(pca, nPC)
  
  ###########################
  #estimate correlation variables vs PCs
  estimate_PCA_corr<-function(pca, pdata, nPC){
    
    #pdata$Basename=make.names(pdata$Basename)
    
    tab<-data.frame(pca$x[,1:nPC]) %>% 
      tibble::rownames_to_column("Basename") %>% 
      tidyr::gather("PC","contrib",-Basename) %>% 
      merge(pdata, by="Basename")
    
    nb_lvl<-tab %>% purrr::map(levels) %>% sapply(length)
    is_num<-tab %>% dplyr::select(-contrib) %>% purrr::map(is.numeric) %>% unlist()
      
    #Wilcox_test for categorical variables with 2 levels
    wilcox<-sapply( names(which(nb_lvl==2) ), function(variable){
        tab %>% group_by(PC) %>% rstatix::wilcox_test(reformulate(variable,"contrib")) %>% pull(p)
    })
    #Kruskal test for categorical variables with more than 2 levels
    kruskal<-sapply( names(which(nb_lvl>2) ), function(variable){
      tab %>% group_by(PC) %>% rstatix::kruskal_test(reformulate(variable,"contrib")) %>% pull(p)
    })
    #cor.test (kendall) for numeric variables
    kendall<-sapply( names(which(is_num==TRUE)), function(variable){
      tab %>% group_by(PC) %>% rstatix::cor_test(!!variable,contrib, method="kendall") %>% pull(p)
    })
    
    dt_pval=cbind(wilcox,kruskal,kendall) %>% t()
    colnames(dt_pval)=paste0( "PC", 1:nPC )
    return(dt_pval)
  }
  
  df_pval<-estimate_PCA_corr(pca, pdata, nPC)
  df_qval<-df_pval %>% as.vector() %>% p.adjust(method='fdr') %>% 
    matrix(ncol=ncol(df_pval), dimnames = list(rownames(df_pval),colnames(df_pval)))
  
  return( list(pvalue=df_pval, qvalue=df_qval, contrib=p, pca=pca))
}
