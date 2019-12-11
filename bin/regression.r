#####################################################################################
#
# Title  : models.r
# Author : CahaisV@iarc.fr
# Package: methylkey
# Date   : 15/11/2019
#
#####################################################################################

#######################
# least square method, return the whole table
m_regression<-function(mval,pdata,variables,method="ls", niter=50, ncore=4){

    regression=NULL
    if(method=="ls"){ 
        regression<-m_leastsquare(mval,pdata,variables)
    }
    if(method=="robust"){ 
        regression<-m_robust(mval,pdata,variables, niter=niter)
    }
    if(method=="logistic"){ 
        regression<-m_logistic(mval,pdata,variables, ncore=ncore)
    }
    return(regression)
}

#######################
# least square method, return the whole table
m_leastsquare<-function(mval,pdata,variables){

    formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
    design<-model.matrix(formula1,data=pdata)
    colnames(design)<-make.names(colnames(design))
    cmtx <- makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )

    fit<-lmFit(mval,design, pdata, ndups=1, method='ls')
    rownames(cmtx)<-colnames(fit)
    fitContrasts=contrasts.fit(fit,cmtx)
    eb=eBayes(fitContrasts)
    chisq <- qchisq(1-eb$p.value,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    table<-topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P")
    rownames(table)<-rownames(mval)

	return(list(table=table, lambda=lambda, pvals=eb$p.value ) )
}

#######################
# least square rosbust method, return the whole table
m_robust<-function(mval,pdata,variables,niter=50){


    formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
    design<-model.matrix(formula1,data=pdata)
    colnames(design)<-make.names(colnames(design))
    cmtx <- makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )

    fit<-lmFit(mval,design, pdata, ndups=1, method='robust', maxiter=niter)
    rownames(cmtx)<-colnames(fit)
    fitContrasts=contrasts.fit(fit,cmtx)
    eb=eBayes(fitContrasts)
    chisq <- qchisq(1-eb$p.value,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    table<-topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P")
    rownames(table)<-rownames(mval)
    
	return(list(table=table, lambda=lambda, pvals=eb$p.value) )
}

#######################
# conditionnal logistic regression, return direcly a toptable
m_logistic<-function(mval,pdata,variables,fdr, ncore=4){

	require(survival)

	tmpval<-merge( pdata[,c("samples",variables) ], t(mval), by.x="samples", by.y="row.names" )
	bcol<-length(c("samples",variables))+1

    #mock result when clogit failed
    mock<-data.frame(matrix( c(0,0,0,0,0,1), nrow=1))
    colnames(mock)<-c("coef","exp(coef)","se(coef)","robust se","z","Pr(>|z|)")

    #clogit parallelized whith mcapply
    logreg<-mclapply( c(bcol:ncol(tmpval)), function(x) {
        tryCatch({
		    as.data.frame( coef(summary( clogit( eval(parse( text=makeFormula(x, variables) )), id=pdata$samples ) )) )
        }, error=function(e){ mock }) 
	}, mc.preschedule=TRUE, mc.cores=ncore )
    logreg<-rbindlist(logreg,fill=TRUE)

	rownames(logreg)<-colnames(tmpval[bcol:ncol(tmpval)])
	table<-data.frame( Odds=logreg$exp,
		        AveExpr=log2( colMeans( m2beta(tmpval[,-c(1:(bcol-1))]) ) ), 
		        t=logreg$z, 
		        P.Value=logreg$Pr, 
		        adj.P.Val=p.adjust(logreg$Pr ,"fdr"), 
		        zvalue=logreg$z )
    table<-table[ !is.na(table$adj.P.Val), ]

    chisq <- qchisq(1-table$P.Value,1)
    lambda<- median(chisq)/qchisq(0.5,1)

    return( list(table=table, lambda=lambda, pvals=table$P.Value ) )
}


#######################
# makeFormula
makeFormula<-function(x, variables){
	
	formule="y~x"
	if( length(variables) > 1){ 
		formule=paste0( "as.numeric(tmpval[,'", variables[1], "']) ~ tmpval[,",x,"] + strata(" , paste( paste0("tmpval[,'",variables[-1], "']"), collapse="+" ), ")")
	}
	if( length(variables) ==1 ) {
		formule=paste0( "as.numeric(tmpval[,'", variables[1], "']) ~ tmpval[,", x, "]" )
	}  
	return(formule)
}