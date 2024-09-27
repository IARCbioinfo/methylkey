

#' m_regression
#' 
#' Display betas in a violin plot by group
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model model
#' @param method ls,ls_ruv,robust,logistic
#' 
#' @return regression
#' 
#' @export
#' 
m_regression<-function(mval,pdata,model,method="ls", niter=50, ncore=4){

    regression=NULL
    if(method=="ls"){ 
        regression<-m_leastsquare2(mval,pdata,model)
    }
    if(method=="ls_ruv"){ 
      regression<-m_leastsquare_ruv(mval,pdata,model)
    }
    if(method=="robust"){ 
        regression<-m_robust(mval,pdata,model, niter=niter)
    }
    if(method=="logistic"){ 
        regression<-m_logistic(mval,pdata,model, ncore=ncore)
    }
    return(regression)
}


#' m_leastsquare
#' 
#' least square method, return the whole table
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model model
#' 
#' @return regression
#' 
m_leastsquare<-function(mval,pdata,model){

    formula1 <- as.formula(model)
    design<-model.matrix(formula1,data=pdata)
    colnames(design)<-make.names(colnames(design))
    cmtx <- limma::makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )
    
    fit<-limma::lmFit(mval,design, pdata, ndups=1, method='ls')
    rownames(cmtx)<-colnames(fit)
    fitContrasts=limma::contrasts.fit(fit,cmtx)
    eb=limma::eBayes(fitContrasts)
    chisq <- qchisq(1-eb$p.value,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    table<-limma::topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P")

    #MetaFile<-data.frame(probeID=rownames(eb$coefficients),Coefficient=eb$coefficients[,1],Stdev=(sqrt(eb$s2.post) * eb$stdev.unscaled)[,1], PValue=eb$p.value[,1])
    table$Coefficient=eb$coefficients[rownames(table),1]
    table$Stdev=(sqrt(eb$s2.post) * eb$stdev.unscaled)[rownames(table),1]

    #goodness
    sst<-rowSums(mval^2)
    ssr<-sst-fit$df.residual*fit$sigma^2
    rsq<-ssr/sst
    table$goodness=rsq[rownames(table)]
    
    return(list(table=table, lambda=lambda ) )
}


# not exported
topTables<-function(eb,x,rsq){
  
  x_table<-limma::topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P", coef=x)
  x_table$Coefficient=eb$coefficients[rownames(x_table),x]
  x_table$Stdev=(sqrt(eb$s2.post) * eb$stdev.unscaled)[rownames(x_table),x]
  x_table$goodness=rsq[rownames(x_table)]
  
  return(x_table)
}


#' m_leastsquare2
#' 
#' least square method, return the whole table
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model model
#' 
#' @return regression
#' 
m_leastsquare2<-function(mval,pdata,model){
  
  formula1 <- as.formula(model)
  design<-model.matrix(formula1,data=pdata)
  colnames(design)<-make.names(colnames(design))
  cmtx <- limma::makeContrasts( contrasts=colnames(design) , levels=colnames(design) )
  
  # lmFit
  fit<-limma::lmFit(mval,design, pdata, ndups=1, method='ls')
  rownames(cmtx)<-colnames(fit)
  fitContrasts=limma::contrasts.fit(fit,cmtx)
  eb=limma::eBayes(fitContrasts)
  
  # lambda
  chisq <- qchisq(1-eb$p.value,1)
  lambda <- apply(chisq, 2, function(x){ median(x)/qchisq(0.5,1) })
  
  #goodness
  sst<-rowSums(mval^2)
  ssr<-sst-fit$df.residual*fit$sigma^2
  rsq<-ssr/sst
  
  # topTables
  topTables<-lapply(colnames(fit), function(x){ topTables(eb,x,rsq) })
  names(topTables)<-colnames(fit)
  
  return(topTables)
}


#' m_leastsquare_ruv
#' 
#' least square method with ruv correction, return the whole table
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model model
#' 
#' @return regression
#' 
#' 
m_leastsquare_ruv<-function(mval,pdata,model){
  
  formula1 <- as.formula(model)
  design<-model.matrix(formula1,data=pdata)
  colnames(design)<-make.names(colnames(design))
  cmtx <- limma::makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )
  fit<-limma::lmFit(mval,design, pdata, ndups=1, method='ls')
  rownames(cmtx)<-colnames(fit)
  fitContrasts=limma::contrasts.fit(fit,cmtx)
  eb=limma::eBayes(fitContrasts)
  lTop<-limma::topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P")
  # Perform RUV adjustment and fit
  ctl <- rownames(mval) %in% rownames(lTop[lTop$adj.P.Val > 0.5,])
  fit <- RUVfit(Y=mval, X=pdata[,factors], ctl=ctl)
  fit2 <- RUVadj(Y=mval, fit=fit)
  chisq <- qchisq(1-fit2$C$F.p,1)
  lambda<-median(chisq)/qchisq(0.5,1)
  table <- topRUV(fit2, sort.by = "F.p", n=nrow(fit2$C))
  return(list(table=table, lambda=lambda, pvals=eb$p.value ) )
}

#' m_robust
#' 
#' least square rosbust method, return the whole table
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model model
#' @param niter number of iteration (50)
#' 
#' @return regression
#' 
m_robust<-function(mval,pdata,model,niter=50){

    formula1 <- as.formula(model)
    design<-model.matrix(formula1,data=pdata)
    colnames(design)<-make.names(colnames(design))
    cmtx <- limma::makeContrasts( contrasts=colnames(design)[2] , levels=colnames(design) )

    fit<-lmFit(mval,design, pdata, ndups=1, method='robust', maxiter=niter)
    rownames(cmtx)<-colnames(fit)
    fitContrasts=limma::contrasts.fit(fit,cmtx)
    eb=limma::eBayes(fitContrasts)
    chisq <- qchisq(1-eb$p.value,1)
    lambda<-median(chisq)/qchisq(0.5,1)
    table<-limma::topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P")
    
    table$Coefficient=eb$coefficients[rownames(table),1]
    table$Stdev=(sqrt(eb$s2.post) * eb$stdev.unscaled)[rownames(table),1]
    
    #goodness
    sst<-rowSums(mval^2)
    ssr<-sst-fit$df.residual*fit$sigma^2
    rsq<-ssr/sst
    
    return(list(table=table, lambda=lambda, pvals=eb$p.value, goodness=rsq) )
}


#' m_logistic
#' 
#' conditionnal logistic regression, return direcly a toptable
#' 
#' @param mval mvalues array
#' @param pdata sampleSheet
#' @param model variables (vector)
#' @param fdr fdr
#' @param ncore number of core (4)
#' 
#' @return regression
#' 
#' @export
#' 
m_logistic<-function(mval,pdata,variables,fdr, ncore=4){

    ##TOREMOVE
    colnames(pdata)[1]<-"samples"
    ##
    
    colnames(mval)<-pdata$samples
    tmpval<-merge( pdata[,c("samples",variables) ], t(mval), by.x="samples", by.y="row.names" )
    bcol<-length(c("samples",variables))+1

    #mock result when clogit failed
    mock<-data.frame(matrix( c(0,0,0,0,0,1), nrow=1))
    colnames(mock)<-c("coef","exp(coef)","se(coef)","robust se","z","Pr(>|z|)")

    #clogit parallelized whith mcapply
    logreg<-mclapply( c(bcol:ncol(tmpval)), function(x) {
        tryCatch({
		    as.data.frame( coef(summary( clogit( eval(parse( text=makeFormula(x, variables) )), id=pdata$samples, robust=FALSE ) )) )
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

    return( list(table=table, lambda=lambda, pvals=table$P.Value, logreg=logreg ) )
}


#' makeFormula
#' 
#' make formula for logistic regression
#' 
#' @param x x
#' @param model variables (vector)
#' 
#' @return regression
#'
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

