#' methyldiff
#'
#' differential analysis,
#' 
#' @param betas matrix of betas
#' @param pdata sample Sheet (dataframe)
#' @param model model to apply with sva (string) eg: "~group+gender"
#' @param case case group
#' @param control control group
#' @param method <limma> ls or robust
#' @param niter <limma> number of iteration for robust method
#' @param ncore <limma> number of core to use for limma
#' @param qval cutoff for significant probes
#' @param sva run sva correction before regression
#' @param modelSVA model to use for batch correction
#'
#' @return toptable for all probes
#'
#' @export
#' 
methyldiff<-function(se=NULL,model="~group",case=NULL,control=NULL,method="ls",niter=50,ncore=2,qval=0.05){
  
  withProgress(message = 'methyldiff', value = 0.1, {
  
    mval=getMvals(se)
    print("mval=")
    print(dim(mval))
    betas=betas_masked_no_na(se)
    pdata=colData(se)
    
    #get manifest
    incProgress(1/8, detail = "get manifest")
    manifest=getAnnotedManifest(metadata(se)$plateform)
    
    save(mval, betas, pdata, method, niter, ncore, qval, manifest, file="/luca/git/shinydmr/test/debug.rda")
    
    #relevel according to case and control
    incProgress(1/8, detail = "relevel")
    grp_g<-strsplit(model,"~|\\+")[[1]][2]
    print(grp_g)
    pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), case)
    pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), control)
    
    #pull first group and sample names into vectors
    grp_g<- pdata %>% dplyr::pull(grp_g) %>% as.factor()
    
    # remove NA in betas again !!
    betas <- betas[rownames(mval),]
    betas <- replaceByMean( betas, groups = grp_g )

    samples <- pdata$samples %>% as.factor()
    #regression analysis
    incProgress(1/8, detail = "regression analysis")
    regression<-m_regression(mval,pdata, model, method=method, niter=niter, ncore=ncore )
    incProgress(1/8, detail = "calculate delta betas")
    deltaBetas<-getDeltaBetas(betas[rownames(regression$table),],grp_g,case,control)*100
    incProgress(1/8, detail = "summarize results")
    table<-data.frame(probeID=rownames(regression$table), regression$table, deltabetas=deltaBetas) %>%
      dplyr::select("probeID","P.Value","adj.P.Val","t","Coefficient","Stdev","deltabetas") %>%
      #tidyr::separate(probeID, sep="_", into=c("probeID","x"),remove=TRUE, fill="right") %>%
      dplyr::mutate(status=ifelse(deltabetas>0,"hyper","hypo")) %>%
      dplyr::mutate(is.sig=adj.P.Val<qval) %>%
      merge(manifest, by="probeID") %>%
      #tidyr::unite("probeID", c(probeID,x),sep="_", na.rm=TRUE, remove=TRUE) %>%
      dplyr::arrange(desc(adj.P.Val)) 
    incProgress(1/8, detail = "done")
    
    #mrs<-MethylResultSet(reac$se,dmps=reac$dmps)
    mrs<- table %>% DataFrame()
    metadata(mrs) = metadata(se)
    mrs@metadata$genome=getGenome(metadata(se)$plateform)
    mrs@metadata$model=model
    mrs@metadata$model.case=case
    mrs@metadata$model.control=control
    mrs@metadata$model.method=method
    mrs@metadata$model.niter=niter
    mrs@metadata$model.ncore=ncore
    mrs@metadata$model.lambda=regression$lambda
    mrs@metadata$model.is.sig=sum(table$is.sig)
    mrs@metadata$model.id=substring(uuid::UUIDgenerate(use.time = FALSE),1,8)
    mrs@metadata$model.date=date=format(Sys.Date(), "%Y%M%d")
    #mrs@metadata$dmrs.qval=qval
    #mrs@metadata$dmrs.pcutoff=pcutoff
    #mrs@metadata$dmrs.maxgap=maxgap
    
    return(mrs)
    
  })
}


#' Batch Correction with SVA
#'
#' This function will do batch correction on mvalues,
#' 
#'
#' @param mval matrix of mvalues
#' @param pdata sample Sheet (dataframe)
#' @param model model to apply with sva (string) eg: "~group+gender"
#' 
#' @return A matrix of batch corrected mvalues
#' 
#' @export
#' 
bc_sva<-function(mval,pdata,model){
  
  print("save")
  save(mval, pdata, model, file="debug.rda")
  print("done")
  
  formula1 <- as.formula(model)
  design<-model.matrix(formula1,data=pdata)
  if(nrow(design) > ncol(mval)){ stop("Missing samples in betas ! ") }
  if(nrow(design) < ncol(mval)){ stop("Missing samples in pdata ! ") }
  
  if (!grepl("\\+",model)){
    sva_m<-sva::sva(mval,design)
  }else{
    formula0 <- as.formula( gsub("~[^+]*\\+","~1+",model) )
    model0<-model.matrix(formula0, data=pdata)
    sva_m<-sva::sva(mval,design,model0)
  }
  
  n.sv=sva_m$n.sv
  message(paste0("number of surrogate variables", sva_m$n.sv))
  if(!sva_m$n.sv) {
    message( "0 surrogate variables have been found, batchcorrection is useless !" )
  }else{
    mval = t(residuals(lm(t(mval)~sva_m$sv)))
  }
  
  return(mval)
}




#' getDeltaBetas
#' 
#' calculate Delta betas
#' 
#' @param betas array of betas values
#' @param group group
#' @param case case
#' @param control control
#' 
#' @return vector
#' 
#' @export
#' 
getDeltaBetas<-function(betas,group, case="TT", control="NT"){
  betas<-betas[,group %in% c(case,control)]
  group<-group[group %in% c(case,control)]
  deltaBetas <- rowMeans(betas[,group==case]) - rowMeans(betas[,group==control])
  return(deltaBetas)
}

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
        regression<-m_leastsquare(mval,pdata,model)
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
    table$goodness=rsq
    
    return(list(table=table, lambda=lambda ) )
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

