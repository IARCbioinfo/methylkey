#####################################################################################
#
# Title  : pca.r
# Author : CahaisV@iarc.fr
# Date   : 28/06/2018
#
#####################################################################################

#########################
##### PCA

### Estimate correlation covs vs. PCA axes

makepca <- function( betas, pdata, out=getwd(), variables=colnames(pdata),  nPC=ncol(betas), id="_raw" ){

	require(data.table)
	require(reshape2)
	require(xtable)
	require(ggplot2)

	invtbetas = 1/t(betas) # n x p required for prcomp
	#mval<-t(beta2m(as.matrix(betas)))
	invtbetas[!is.finite(invtbetas)]<-min(invtbetas[is.finite(invtbetas)])

	#variance for a probe should be > 0
	sel<-which(apply(invtbetas, 2, var)==0)

	if (length(sel)) { 
		pca = prcomp(invtbetas[,-sel],scale=T) 
	} else { 
		pca = prcomp(invtbetas,scale=T) 
	}

	# plot proportion of var explained
	propVarExpl = summary(pca)$importance[2,]

        ###########################
	#plot contributions
	plot_PCA_contribution <- function( pca, nPC = 10, title = "" ){
		v_var <- summary( pca )$importance[ 2, 1:nPC ] 
		return( ggplot(data.table( PCA_axis = names( v_var ), val = v_var )) 
		+ geom_col(aes( x=reorder( PCA_axis, 1:nPC ), y=val )) 
		+ labs( title = title, subtitle = "PCA contribution", x = "PC", y = "Prption expl. var." )
		+ theme( plot.title = element_text(lineheight=.8, face="bold") ) 
		)
	}

	p<-plot_PCA_contribution(pca, nPC)
	jpeg(paste0(out, "/pca_contributions_", id, ".jpg"))
	plot(p)
	dev.off()

	###########################
	#estimate correlation variables vs PCs
	estimate_PCA_corr<-function(pca, pdata, nPC, variables){
		
		tab<-data.table( cbind( pdata, as.data.table( pca$x[,1:nPC] ) ))
		PCs = paste0( "PC", 1:nPC )
		dt_pval <- data.table( PCA_dim = unlist( lapply( PCs, rep, length( variables ) ) ), Variables = rep( variables, length( PCs ) ), p_val= rep( 2, length(PCs ) ) )

		for (PC in PCs){
			for ( variable in variables ){
				if ( tab[ , is.factor( get(variable) ) ]  ){
					if ( length( levels( tab[ ,  get(variable)  ] ) ) == 2 ){
						pv <- tab[ , wilcox.test(  get(PC) ~ get(variable)  ) ]$p.value
					}else if ( length( levels( tab[ ,  get(variable)  ] ) ) > 2 ) {
						pv <- tab[ , kruskal.test(  get(PC) ~ get(variable)  ) ]$p.value
					}else { pv<-1 }

				}else if ( tab[ , is.numeric( get(variable) ) ] ){
					pv <- tab[ , cor.test(  get(PC), get(variable) , method = "kendall" ) ]$p.value
				}
				dt_pval[ PCA_dim == PC & Variables == variable, ]$p_val <- pv
			}
		}

		dt_pval[ , adj_pval := p.adjust( p_val, method = "fdr" ) ]

	  	dt_pval[ , class_pval := "NS" ]
	  	dt_pval[ adj_pval < 0.05, class_pval := "." ]
	  	dt_pval[ adj_pval < 0.01, class_pval := "*" ]
	  	dt_pval[ adj_pval < 0.001, class_pval := "**" ]
	  	dt_pval[ adj_pval < 0.0001, class_pval := "***" ]
	  	dt_pval[ adj_pval < 0.00001, class_pval := "****" ]
		dt_pval[ , PCA_dim2 := as.numeric( gsub( "PC", "", PCA_dim ) ) ]

		return(data.table(dt_pval))
	}

	dt_pval<-estimate_PCA_corr(pca, pdata, nPC, variables)
	
	mt_pval<-acast(dt_pval, Variables~PCA_dim, value.var="p_val")
	html_pval<-print( xtable(mt_pval), type="html", print.results=FALSE)
	
	mt_qval<-acast(dt_pval, Variables~PCA_dim, value.var="adj_pval")
	html_qval<-print( xtable(mt_qval), type="html", print.results=FALSE)
	
	return( list(pvalue=html_pval, qvalue=html_qval))
}

