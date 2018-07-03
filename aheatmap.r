library(data.table)
bpa<-fread("bpa.txt")
bpf<-fread("bpf.txt")
bps_h<-fread("bps_h.txt")

foo<-merge(bpa[,c(4,12)], bpf[,c(4,12)], by="cpg", all=T)
colnames(foo)<-c("cpg","bpa","bpf")
foo<-foo[ is.na(bpa) | is.na(bpf), ]

foo<-merge(foo, bps_h[,c(4,12)], by="cpg", all=T)
colnames(foo)<-c("cpg","bpa","bpf","bps_h")
foo<-foo[ (is.na(bpa) | is.na(bps_h)) & (is.na(bpf) | is.na(bps_h)),  ]

foo[ , max:=pmax(abs(bpa),abs(bpf),abs(bps_h), na.rm=T) ]
#fwrite(foo,file="bpafs.txt")


library(NMF)
betas<-read.table("betas",header=T, sep="\t")
pdata<-fread("pdata")
foo<-fread("heatmap_BPAFS")
foo[ , max:=pmax(abs(bpa_deltaBetas),abs(bpf_deltaBetas),abs(bps_deltaBetas), na.rm=T) ]

sel<-foo[ order(max, decreasing=T) , cpgs ]
jpeg("heatmap_BPAFS.jpeg", width=2800, height=2000, res = 200) 
aheatmap(betas[sel,], annCol=pdata$Sample_condition, annRow=foo[sel]$group, labCol=pdata$Sample_ID)
dev.off()



