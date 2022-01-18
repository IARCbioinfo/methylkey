#' readSampleSheet
#' 
#' Display betas in a violin plot by group
#' 
#' @param path file location
#' @param sample sample column name
#' @param barcode barcode column name
#' @param groups groups columns names (vector)
#' @param sep field separator ('\\t')
#' 
#' @return dataframe
#' 
#' @export
#' 
readSampleSheet<-function(path, samples=NULL, barcode=NULL, groups=NULL, sep="\t"){

    #1- read file
    pdata<-readr::read_delim(path,delim=sep,trim_ws=TRUE)

    if( ncol(pdata)<2 ){ stop( paste0("Oups ! Check if your separator in pdata is really '", sep, "' in your sample sheet.") ) }
    
    #identify columns indexes
    samples_<-which( tolower( colnames(pdata) ) %in% c("samples","sample","samples_names","sample_names","sample_name","sample_id","id","ids") )
    barcode_<-which( tolower( colnames(pdata) ) %in% c("barcode","basename") )
    sentrix_id_<-which( tolower( colnames(pdata) ) %in% c("sentrix_id") )
    sentrix_position_<-which( tolower( colnames(pdata) ) %in% c("sentrix_position") )

    #2 samples column
    if( !is.null(samples) ){ 
        colnames(pdata)[ which( colnames(pdata)==samples ) ]<-"samples"
    }else if (is.numeric(samples)){
        colnames(pdata)[ samples ]<-"samples"
    }else if ( length(samples_) ) {
        colnames(pdata)[samples_]="samples"
    }else { stop( "Oups, sorry, cannot determine samples names column. Try to use --samples <column_name>" )   }
    if ( !base::all( !duplicated(pdata$samples) ) ) { stop("Oups ! there is duplicated samples name in your sample sheet.") }
    message( paste0("reading samples sheet : ", nrow(pdata), " samples") )
    #3- Basename columns
    if( !is.null(barcode) ){ 
        colnames(pdata)[ which( colnames(pdata)==barcode ) ]<-"Basename"
    }else if (is.numeric(barcode)){
        colnames(pdata)[ barcode ]<-"Basename"
    }else if ( length(barcode_) ) {
       colnames(pdata)[barcode_]="Basename"
    }else if ( length(sentrix_id_) == 1 & length(sentrix_position_) ==1 ){
      if ( nchar( as.character( pdata[sentrix_position_][1,1] ) ) == 3 ) { 
            pdata[sentrix_position_]<-paste0(unlist(pdata[sentrix_position_]),"C01")
      }
      pdata$Basename<-paste0( pdata[,sentrix_id_], "_" , pdata[,sentrix_position_] )
    }else { stop( "Oups, sorry, cannot determine barcode column. Try to use --barcode <column_name>" )}

    #4- convert discrete variables to factor
    pdata <- pdata %>% mutate_if(~ length(unique(.)) <= length(.)/8 , as.factor )
      
    #5- check groups
    for (group in groups) {
        if ( ! group %in% colnames(pdata) & !is.numeric(group) ){
            stop(paste0("Oups ! Check if '", group, "' is a column in your sample sheet"))
        }
    }

    message( paste0("reading samples sheet : done ...") )
    return(pdata)
}
