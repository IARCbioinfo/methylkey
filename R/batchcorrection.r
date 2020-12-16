bc_sva<-function(mval,pdata,model){

  library(sva)
  
  formula1 <- as.formula(params$model)
  design<-model.matrix(formula1,data=pdata)
  if(nrow(design) > ncol(mval)){ stop("Missing samples in betas ! ") }
  if(nrow(design) < ncol(mval)){ stop("Missing samples in pdata ! ") }
  
  if (!grepl("\\+",model)){
    sva<-sva(mval,design)
  }else{
    #formula0 <- as.formula(paste("~ 1+", paste(variables[-1],collapse="+")))
    formula0 <- gsub("~[^+]*\\+","~1+",formula1)
    model0<-model.matrix(formula0, data=pdata)
    sva<-sva(mval,design,model0)
  }

  n.sv=sva$n.sv
  message(paste0("number of surrogate variables", sva$n.sv))
  if(!sva$n.sv) {
    message( "0 surrogate variables have been found, batchcorrection is useless !" )
  }else{
    mval = t(residuals(lm(t(mval)~sva$sv)))
  }

  return(mval)
}