#####################################################################################
#
# Title  : missingValues.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 20/06/2018
#
#####################################################################################


#select probes with percentage of missing values > CpGlimit
CpGNAexcl<-function(betas,nalimit=0.2){
	
	print(paste0("remove probes with number of missing values > ", nalimit ))

	NArow<-apply(betas,1,function(i) sum(is.na(i)))
	exclprob<-rownames(betas)[NArow>nalimit*ncol(betas)]
	return(exclprob)
}


#impute missing values
imputeNA<-function(betas,nalimit){

	print("impute missing values with pamr")

	require(pamr)
	#impute missing values
	foo<-list( x=as.matrix(betas), y=colnames(betas) )
	betas<-pamr.knnimpute(foo, rowmax=nalimit , colmax=0.95)$x
	return(betas)
}


#replace missing values by mean for each group
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
