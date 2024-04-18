#' Define a Class Union
#'
#' This function allows you to define a class union, which is a special type of class
#' that includes objects from multiple classes. A class union is useful when you want
#' to create a class that can accept objects of different classes as its instances.
#'
#' @keywords classes
#' 
#' @rdname setClassUnion
setClassUnion("List_OR_NULL", c("list", "NULL"))

#' MethylResultSet Class
#'
#' A class for storing and working with results from DNA methylation analysis.
#' The class contains the following slots:
#' \describe{
#'   \item{dmps}{A data frame (DataFrame) containing differential methylation positions (DMPs).}
#'   \item{dmrs}{A list (List_OR_NULL) containing differential methylation regions (DMRs).}
#'   \item{pca}{A list containing principal component analysis (PCA) results.}
#'   \item{metadata}{A list containing metadata or additional information about the results.}
#' }
#'
#' @docType class
#' @name MethylResultSet
#' @aliases MethylResultSet-class
#' @author Your Name
#'
#' @slot dmps A data frame (DataFrame) containing differential methylation positions (DMPs).
#' @slot dmrs A list (List_OR_NULL) containing differential methylation regions (DMRs).
#' @slot pca A list containing principal component analysis (PCA) results.
#' @slot metadata A list containing metadata or additional information about the results.
#'
#' @rdname MethylResultSet
#'
#' @export
MethylResultSet<-setClass(
  
  "MethylResultSet",
  representation(
    manifest="DFrame",
    dmps="List_OR_NULL",
    lambda="numeric",
    metadata="list"
  ),
  prototype(
    manifest=S4Vectors::DataFrame(),
    dmps=list(),
    lambda=c(),
    metadata=list()
  )
  
)

# not exported
topTables<-function(eb,x,rsq){
  
  x_table<-limma::topTable(eb, adjust="BH", number=Inf, p=1, sort.by="P", coef=x)
  x_table$Coefficient=eb$coefficients[rownames(x_table),x]
  x_table$Stdev=(sqrt(eb$s2.post) * eb$stdev.unscaled)[rownames(x_table),x]
  x_table$goodness=rsq[rownames(x_table)]
  
  return(x_table)
}

#' Create a MethylResultSet Object
#'
#' Create an instance of the MethylResultSet class to store DNA methylation analysis results.
#'
#' @param se A MethylSet object containing DNA methylation data.
#' @param dmps A data frame containing differential methylation positions (DMPs).
#' @param dmrs A list of data frames containing differential methylation regions (DMRs).
#' @param pca A list containing principal component analysis (PCA) results.
#'
#' @return An object of class MethylResultSet.
#'
#' @seealso \code{\link{MethylResultSet-class}}
#'
#'
#' @rdname MethylResultSet
#'
#' @export
MethylResultSet <- function(se,model,intercept, method="ls")
{
  
  manifest = methylkey::getManifest( metadata(se)$plateform )
  pdata = data.frame(colData(se))
  mval = getMvals(se)
  model = tolower(model)
  
  grp_g<-strsplit(model,"~|\\+")[[1]][2]
  pdata[,grp_g] <- relevel(as.factor(unlist(pdata[,grp_g])), intercept)
  
  formula1 <- as.formula(model)
  design<-model.matrix(formula1,data=pdata)
  colnames(design)<-make.names(colnames(design))
  cmtx <- limma::makeContrasts( contrasts=colnames(design) , levels=colnames(design) )
  
  # lmFit
  fit<-limma::lmFit(mval,design, pdata, ndups=1, method=method)
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
  
  dmps<-lapply(colnames(fit)[-1], function(x){ topTables(eb,x,rsq) })
  names(dmps)<-colnames(fit)[-1]
  
  ### J'en suis ici, ajouter les delta betas !
  
  for( x in names(dmps)){
    contrast<-paste0(gsub(grp_g,"",x),"_vs_",intercept)
    add_deltabetas = contrast %in% names(rowData(methM))
    dmps[[x]] <- dmps[[x]] %>%
      { if(add_deltabetas) dplyr::mutate(., deltabetas = rowData(methM)[ rownames(dmps[[x]]), contrast ] ) else .  } %>% # add delta betas
      { if(add_deltabetas) dplyr::mutate(., status=ifelse(deltabetas>0,"hyper","hypo")) else .  } %>%
      tibble::rownames_to_column("Probe_ID") %>%
      dplyr::select(-logFC) %>%
      dplyr::arrange(Probe_ID)
  }
  
  metadata=metadata(se)
  metadata$model=model
  metadata$intercept=intercept
  
  index_order=match(dmps[[1]]$Probe_ID,manifest$Probe_ID,nomatch = 0)
  manifest=manifest[index_order,]
  
  new("MethylResultSet",
      manifest=S4Vectors::DataFrame( manifest ),
      dmps=dmps,
      lambda=lambda,
      metadata=metadata
  )
  
}

#' Create generic function for retrieving differentially methylated positions (DMPs)
#'
#' This function sets up a generic function, `getDMPs`, for retrieving differentially methylated positions (DMPs) from a MethylResultSet object.
#'
#' @title Create generic function for retrieving DMPs
#' @param x A MethylResultSet object.
#' @param group The group number for which to retrieve DMPs (default is 1).
#' @export
setGeneric("getDMPs", function(x,group=1) 
  standardGeneric("getDMPs") )

setMethod("getDMPs", "MethylResultSet",
  function(x,group=1){
    return(dplyr::left_join(x@dmps[[group]],as_tibble(x@manifest),by="Probe_ID"))
})

#' Create generic function for retrieving differentially methylated position (DMP) ranges
#'
#' This function sets up a generic function, `getDMPranges`, for retrieving ranges of differentially methylated positions (DMPs) from a MethylResultSet object.
#'
#' @title Create generic function for retrieving DMP ranges
#' @param x A MethylResultSet object.
#' @param group The group number for which to retrieve DMPs (default is 1).
#' @param q The threshold for adjusted p-value (default is 0.5).
#' @export
setGeneric("getDMPranges", function(x,group=1,q=0.5) standardGeneric("getDMPranges") )

setMethod("getDMPranges", "MethylResultSet", function(x,group=1,q=0.05){
  df = left_join(as.data.frame(x@manifest), x@dmps[[group]], by="Probe_ID") %>%
    dplyr::filter(adj.P.Val < q)
  GenomicRanges::makeGRangesFromDataFrame(df, start.field="pos", end.field="pos", keep.extra.columns = TRUE )
})

#' Create generic function for retrieving methylation analysis results
#'
#' This function sets up a generic function, `getResults`, for retrieving methylation analysis results from a MethylResultSet object.
#'
#' @title Create generic function for retrieving methylation analysis results
#' @param x A MethylResultSet object.
#' @param index Index or indices specifying the samples for which to retrieve results.
#' @param tools Character vector specifying the analysis tools to use (default is c("dmrcate", "ipdmr")).
#' @param ... Additional arguments to be passed to the analysis functions.
#' @export
setGeneric("getResults", function(x,index,tools=c("dmrcate","ipdmr"),...) 
  standardGeneric("getResults") )

setMethod("getResults", "MethylResultSet",
  function(x,index,tools=c("dmrcate","ipdmr"),...){
    
    result<-getDMPs(mrs,index)
    result<-add_dmrcate(result,...)
    result<-add_ipdmr(result,...)
    return(result)
})            

#' Convert results to long format
#'
#' This function takes a data frame containing results from different tools 
#' (e.g., dmrcate, ipdmr, combp, dmrff) and converts it to long format.
#'
#' @title Convert results to long format
#' @param df A data frame containing the results.
#' @return A data frame in long format.

#' @export
setGeneric("resultsToLong", function(df) standardGeneric("resultsToLong") )

setMethod("resultsToLong", "data.frame", function(df){
            
  df %>% pivot_longer(cols = any_of(c("dmrcate","ipdmr","combp","dmrff")), names_to="dmrtool", values_to="dmrs" ) %>%
    separate(dmrs, into=c("ID","fdr","no.cpgs"), sep=":") %>%
    mutate(fdr=as.numeric(fdr), no.cpgs=as.numeric(no.cpgs)) %>%
    dplyr::filter( !is.na(fdr) )
})

#' Create generic function for converting results to long format
#'
#' This function sets up a generic function, `resultsToLong`, for converting results to long format.
#'
#' @title Create generic function for converting results to long format
#' @param df A data frame containing results.
#' @export
setGeneric("add_dmrcate", function(x,fdr=0.05,pcutoff=0.02,maxgap=1000,genome="hg38") 
  standardGeneric("add_dmrcate") )

setMethod("add_dmrcate", "data.frame",
  function(x,fdr=0.05,pcutoff=0.02,maxgap=1000,genome="hg38"){
    
    searchDMR_dmrcate(x,fdr=fdr,pcutoff=pcutoff,maxgap=maxgap,genome=genome) %>%
      unite("dmrcate",seqnames,start,end,sep = "-") %>%
      unite("dmrcate",dmrcate,HMFDR,no.cpgs,sep=":")       
})

#' Create generic function for adding IPDMR results to a data frame
#'
#' This function sets up a generic function, `add_ipdmr`, for adding IPDMR results to a data frame.
#'
#' @title Create generic function for adding IPDMR results
#' @param x A data frame.
#' @param maxgap Maximum gap size for IPDMR regions (default is 1000).
#' @export
setGeneric("add_ipdmr", function(x,maxgap=1000) 
  standardGeneric("add_ipdmr") )

setMethod("add_ipdmr", "data.frame",
  function(x,maxgap=1000){
    
    dmrs<-searchDMR_ipdmr(x, maxgap=maxgap) %>%
      separate_rows(probe, sep = ";") %>% 
      dplyr::rename(Probe_ID=probe) %>% 
      unite("ipdmr",chr,start,end,sep = "-") %>%
      unite("ipdmr",ipdmr,fdr,nprobe,sep = ":") %>%
      mutate(ipdmr = ifelse(!startsWith(ipdmr, "chr"), paste0("chr", ipdmr), ipdmr)) %>%
      dplyr::select(-p)
    
    left_join(x,dmrs,by="Probe_ID")
})


#' \code{getDMRs} Generic Function
#'
#' Description of the generic function \code{getDMRs}.
#'
#' @title getDMRs Generic Function
#' @param x An object.
#' @return A result.
#' @export
setGeneric("getDMRs", function(x) 
  standardGeneric("getDMRs") )

#' getDMRs Method for MethylResultSet Objects
#'
#' Description of the method for extracting DMRs from a MethylResultSet object.
#'
#' @rdname getDMRs
#' @param x An object of class \code{MethylResultSet}.
#' @return A data frame containing unified DMRs.
#' @export
setMethod("getDMRs", "MethylResultSet",
  function(x){
    
    lapply(names(x@dmrs), function(n) .dmrs.unify( x, x@dmrs[[ n ]], n ) ) %>%
      bind_rows()
  })

#' \code{.dmrs.unify} Function
#'
#' Description of the function \code{.dmrs.unify}.
#'
#' @title .dmrs.unify Function
#' @param x An object.
#' @param dt A data frame containing DMR information.
#' @param t A character indicating the DMR tool.
#' @return A data frame with unified DMR information.
#' @export
.dmrs.unify<-function(x,dt,t){
  if(nrow(dt)==0) return(data.frame())
  dt %>% as_tibble() %>% 
    dplyr::select(probes,ID,chr,start,end,no.cpgs,fdr) %>%
    tidyr::unnest(probes) %>%
    left_join(data.frame(x@dmps) %>% dplyr::select("Probe_ID","pos","deltabetas","status"), by="Probe_ID") %>%
    mutate(dmrtool=t)
}
    
