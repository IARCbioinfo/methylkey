violin_plot1<-function(betas, group, numPositions=1000){
  
  require(ggplot2)
  require(data.table)
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-cbind(group=as.character(group),data.table(t(cpg)) )
  cpg<-melt( cpg, id.vars="group" )
  colnames(cpg)<-c("group","probeID","betas")
  
  #violin plots
  p1<-ggplot(cpg, aes(x=group, y=betas, fill=group)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("")
  
  return(p1)
}


violin_plot2<-function(betas, annot, feature, numPositions=1000){
  
  require(ggplot2)
  require(data.table)
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-melt(cpg)
  colnames(cpg)<-c("probeID","Basename","betas")
  cpg<-merge(annot,cpg,by="probeID")
  
  #violin plots
  p1<-ggplot(cpg, aes_string(x=feature, y="betas", fill=feature)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("")
  
  return(p1)
}


violin_plot3<-function(betas, annot, group, feature, numPositions=1000){
  
  require(ggplot2)
  require(data.table)
  
  cpg<-betas[runif(numPositions,1,nrow(betas)),]
  cpg<-cbind(group=as.character(group),data.table(t(cpg)) )
  cpg<-melt( cpg, id.vars="group" )
  colnames(cpg)<-c("group","probeID","betas")
  cpg<-merge(annot,cpg,by="probeID")
  
  #violin plots
  p1<-ggplot(cpg, aes_string(x=feature, y="betas", fill=feature)) + geom_violin() + 
    theme(axis.text.x=element_text(angle=50,hjust=1, size=14)) + 
    guides(fill =FALSE) +
    ylab("Methylation %") + xlab("") +
    facet_wrap(~group)
  
  return(p1)
}

