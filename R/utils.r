#####################################################################################
#
# Title  : utils.r
# Author : CahaisV@iarc.fr, Novoloacaa@student.iarc.fr
# Package: methylkey
# Date   : 17/01/2018
#
#####################################################################################


#######################
# get Delta Betas

getDeltaBetas<-function(betas,group, case="TT", control="NT"){
  betas<-betas[,group %in% c(case,control)]
  group<-group[group %in% c(case,control)]
  deltaBetas <- rowMeans(betas[,group==case]) - rowMeans(betas[,group==control])
  return(deltaBetas)
}

#getDeltaBetas<-function(betas,group){

#	if( is.null(nrow(betas)) ){
#		means<-by( as.matrix(betas,ncol=1), group, colMeans)
#	} else {
#		means<-by(t(betas), group, colMeans)
#	}

#	delta=list()
#	for( level1 in levels(group) ) {
#		for( level2 in levels(group) ) {
#			if (level1 != level2){
#				delta[[paste0(level1,"vs",level2)]]<-means[[level2]]-means[[level1]]
#			}
#		}
#	}
#	return(do.call(cbind, delta))
#}


#######################
# id maker
idmaker <- function(x)
{
    max.val = x*100
    count <- nchar(as.character(max.val))                       # find out how many 'numbers' each ID will have after the letter
    size <- paste("%0",count,"d",sep="")                        # set the variable to be fed into 'sprintf' to ensure we have leading 0's
    lets <- toupper(sample(letters,x, replace=T))               # randomising the letters 
    nums <- sprintf(size,sample(1:max.val)[1])                  # randominsing the numbers, and ensuing they all have the same number of characters
    ids <- paste(lets,nums,sep="")                              # joining them together
    return(ids)
}

#######################
# formatParameters
autoformat<-function(params){

	#param is numeric
	params<-unlist(strsplit(params,","))
	if ( !is.na( as.numeric( params ) ) ) {
		return( as.numeric(params) )
	}
	return(params)
}



###################
#format chr name to numeric
ChrNameToNumeric<-function(chr){
	
	chr<-gsub("chr","",chr)
	chr<-gsub("X","23",chr)
	chr<-gsub("Y","24",chr)
	return(as.numeric(chr))
}



##################
#appendEnv
appendEnv<-function(e1,e2){
  
  e1name = deparse(substitute(e1))
  e2name = deparse(substitute(e2))
  listE1 = ls(e1)
  listE2 = ls(e2)
  for(v in listE2) {
    if(v %in% listE1) warning(sprintf("Variable %s is in e1, too!", v))
    e1[[v]] = e2[[v]]
  }
  return(e1)
}




