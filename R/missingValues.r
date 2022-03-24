#' CpG with too many NA
#'
#' select probes with percentage of missing values superior to CpGlimit
#' 
#' @param betas matrix of betas
#' @param nalimit maximum proportion of NA accepted
#' 
#' @return List of probes to exclude
#' 
#' @export
#' 
CpGNAexcl<-function(betas,nalimit=0.2){

	NArow<-apply(betas,1,function(i) sum(is.na(i)))
	exclprob<-rownames(betas)[NArow>nalimit*ncol(betas)]
	return(exclprob)
}


#' impute missing values
#'
#' impute missing values with parmr for probes with less than 'nalimit' proportion of NA values
#' 
#' @param betas matrix of betas
#' @param nalimit maximum proportion of NA accepted
#' 
#' @return betas matrix
#' 
#' @export
#' 
imputeNA<-function(betas,nalimit){

	require(pamr)
	foo<-list( x=as.matrix(betas), y=colnames(betas) )
	betas<-pamr::pamr.knnimpute(foo, rowmax=nalimit , colmax=0.95)$x
	return(betas)
}


#' replaceByMean
#'
#' Replace missing values by mean for each group
#' 
#' @param betas matrix of betas
#' @param group group name
#' 
#' @return betas matrix
#' 
#' @export
#' 
replaceByMean<-function(betas,groups){
	
	print("replace missing values by means of groups")
	
	#select rows containing missing values
	probNAcont<-which(apply(betas,1,function(i) sum(is.na(i)))>0)

	#calculate means by group
	meanBetas<-function(betas){
		rmeans<-matrix( rowMeans(betas,na.rm=T), ncol=ncol(betas), nrow=nrow(betas) )
		betas[is.na(betas) ] <- rmeans[is.na(betas)]
		return( betas )
	}

	#replace missing values by means
	for(group in levels(groups)){
		sel<-groups==group
		if(!sum(sel)>1){ stop("You should have at least 2 samples by group !") }
		betas[probNAcont,sel]<-meanBetas(betas[probNAcont,sel])
	}
	return(betas)
}


