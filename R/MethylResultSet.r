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
MethylResultSet<-setClass(
  
  "MethylResultSet",
  representation(
    dmps="DataFrame",
    dmrs="List_OR_NULL",
    pca="list",
    metadata="list"
  ),
  prototype(
    dmps=new("DFrame"),
    dmrs=list(),
    pca=list(),
    metadata=list()
  )
  
)

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
MethylResultSet <- function(se, dmps, dmrs=NULL, pca=NULL)
{
  new("MethylResultSet",
      dmps=DataFrame(dmps),
      dmrs=lapply( as.list(dmrs), DataFrame ),
      pca=as.list(pca),
      metadata=metadata(se)
  )
}

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
  dt %>% as_tibble() %>% 
    dplyr::select(probes,ID,chr,start,end,no.cpgs,fdr) %>%
    tidyr::unnest(probes) %>%
    left_join(data.frame(x@dmps) %>% dplyr::select("probeID","pos","deltabetas","status"), by="probeID") %>%
    mutate(dmrtool=t)
}
    
