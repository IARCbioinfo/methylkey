#! /usr/bin/env nextflow

/*vim: syntax=groovy -*- mode: groovy;-*- */

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help          		= null
params.pdata           		= null
params.input      		= null
params.html	  		= "./"
params.out         		= "report_html"
params.samples			= null
params.groups       		= null
params.barcode                  ="barcode"
params.platform			= ""
params.pipeline   		= "minfi"
params.missing			= "mean"
params.normalize		= "funnorm"
params.nalimit  		= "0.2"
params.genome			= "hg19"
params.filters           	= null
params.regions			= ""
params.violin			=FALSE

log.info ""
log.info "-----------------------------------------------------------------------"
log.info "  Methylkey 1.2 : Use methylkey pipeline for methylation data analysis "
log.info "-----------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "-------------------SOMATIC -----------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/methylkey.nf --pdata pdata.txt --idat idat_rep/ --samples samples_names --groups samples_groups"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--pdata                FILE                 Sample sheet"
    log.info "--input                FOLDER/FILE          Folder containing idat files (or file with beta values if platform=matrix)"
    log.info "--samples              STRING/INTEGER       Column name or column index for samples names"
    log.info "--groups               FILE                 Column name or column index for samples groups (eg : Case/Control)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "--html                 STRING               Name for output report (default=report.html)"
    log.info "--out		     STRING               Name for output directory (default=report_files)"
    log.info "--platform             STRING               IlluminaHumanMethylation450k, IlluminaHumanMethylationEPIC or matrix"
    log.info "--genome		     STRING		  Reference genome (default hg19)
    log.info "--barcode              STRING/INTEGER       Column name or column index for array barcodes (default=barcode)"
    log.info "--pipeline             STRING               minfi or methylumi (default=minfi)"
    log.info "--normalize            STRING      	  Default is funnorm with minfi, bmiq with methylumi, use none to desactivate"
    log.info "--filters		     FILE		  File containing list of probes to filter (one probe id by line)"
    log.info "--nalimit		     DOUBLE               Maximum fraction of missing values for a probe"
    log.info "--missing		     STRING               keep, mean (calculate mean of betas by group) or impute (with pamr)"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Flags:"
    log.info "--violin					  Draw the violin plots of beta values by groups
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Cell Mixture Composition:"
    log.info "--cell houseman                             Houseman, only with minfi package, only for blood cells"
    log.info "--nbc		     INTEGER              Number of cell type (default=6)"
    log.info ""
    log.info ""
    exit 1
} else {


pdata = file(params.pdata)


/* Software information */
log.info ""
log.info "pdata           	= ${params.pdata}"
log.info "input		      	= ${params.input}"
log.info "samples  		= ${params.samples}"
log.info "groups	       	= ${params.groups}"
log.info "barcode        	= ${params.barcode}"
log.info "platform           	= ${params.platform}"
log.info "genome           	= ${params.genome}"
log.info "report                = ${params.html}"
log.info "report files        	= ${params.out}"
log.info "pipeline  		= ${params.pipeline}"
log.info "normalize         	= ${params.normalize}"
log.info "filters           	= ${params.filters}"
log.info "NA cutoff    	  	= ${params.nalimit}"
log.info "missing values 	= ${params.missing}"
log.info "violin plot		= ${params.violin}"
log.info ""
log.info "CMC   		= ${cell}"
log.info "number of cell type	= ${nbc}"
log.info ""
}


process run_strelka {

	input:

	file pdata

     shell:
     '''
     Rscript methylkey.r --pdata !{params.pdata} --idat !{params.idat} --html !{params.html} --out !{params.out} --samples !{params.samples} --groups !{params.groups} --pipeline !{params.pipeline} --nalimit !{params.nalimit} --filters !{params.filters}
     '''
}




