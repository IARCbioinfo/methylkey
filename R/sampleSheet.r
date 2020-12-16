#####################################################################################
#
# Title  : pdata.r
# Author : CahaisV@iarc.fr
# Package: methylkey
# Date   : 11/11/2019
#
#####################################################################################

message("loading sampleSheet ...")
#######################
readSampleSheet<-function(path, samples=NULL, barcode=NULL, groups=NULL, sep="\t"){

    #1- read file
    pdata<-read.table(path, sep=sep, header=T)

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

    pdata$Basename<-as.character(pdata$Basename)

    #4- convert discrete variables to factor
    for (variable in colnames(pdata)){ 
        if (length(levels(factor(pdata[,variable]))) <= nrow(pdata)/8) 
        { 
  	        pdata[,variable]<-as.factor(pdata[,variable]) 
        }  
    }

    #5- check groups
    for (group in groups) {
        if ( ! group %in% colnames(pdata) & !is.numeric(group) ){
            stop(paste0("Oups ! Check if '", group, "' is a column in your sample sheet"))
        }
    }

    message( paste0("reading samples sheet : done ...") )
    return(pdata)
}

removeBadSamples<-function(pdata, path){

    bad<-read.table(path)$V1
    pdata<-pdata[ !pdata[,samples] %in% bad, ]
    message( paste0( "reading samples sheet : ", length(bad), " samples listed in ", path, " have been removed !") )
    return(pdata)
}
