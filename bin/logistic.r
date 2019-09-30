#####################################################################################
#
# Title  : logistic.r
# Author : CahaisV@iarc.fr
# Package: methylkey
# Date   : 18/09/2019
#
#####################################################################################


#######################
# logistic regression, return direcly a toptable
logistic<-function(mval,pdata,samples,variables,fdr){

	pdata$samples<-samples
	tmpval<-merge( pdata[,c("samples",variables) ], t(mval), by.x="samples", by.y="row.names" )
	bcol<-length(c("samples",variables))+1
	
	logreg<-sapply( c(bcol:ncol(tmpval)), function(x) {
		coef(summary( glm( eval(parse( text=makeFormula(x, variables) )), family='binomial' ) ))
	}[2,] )

	rownames(logreg)<-c("Estimate","StdError","zvalue", "Pr")
	colnames(logreg)<-colnames(tmpval[bcol:ncol(tmpval)])
	logreg<-as.data.frame(t(logreg))
	toptable<-data.frame( Estimate=logreg$Estimate,
		        AveExpr=log2( colMeans( m2beta(tmpval[,-c(1:(bcol-1))]) ) ), 
		        t=logreg$zvalue, 
		        P.Value=logreg$Pr, 
		        adj.P.Val=p.adjust(logreg$Pr ,"fdr"), 
		        zvalue=logreg$zvalue )
	  
	return(toptable)
}

#######################
# makeFormula
makeFormula<-function(x, variables){
	
	formule="y~x"
	if( length(variables) > 1){ 
		formule=paste0( "tmpval[,2] ~ tmpval[,",x,"] + " , paste( paste0("tmpval[,'",variables[-1], "']"), collapse="+" ))
	}
	if( length(variables) ==1 ) {
		formule=paste0( "tmpval[,2] ~ tmpval[,", x, "]" )
	}  
	return(formule)
}



