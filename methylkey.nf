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
params.tool                     = "methylkey"
params.pdata           		= null
params.idat      		= null
params.out         		= "methylkey"
params.samples			= null
params.groups       		= null
params.barcode                  = "barcode"
params.platform			= ""
params.pipeline   		= "minfi"
params.missing			= "mean"
params.normalize		= "funnorm"
params.nalimit  		= "0.2"
params.genome			= "hg19"
params.filters           	= null
params.regions			= ""
params.violin			= false
params.nbc                      = null
params.cell                     = 6

params.meth                     = "methylkey/meth.rdata"
params.correction               = "sva"
params.batch                    = "no"
params.what                     = "betas"
params.variables                = null

params.model                    = "None"
params.case                     = null
params.control                  = null
params.fdr                      = 0.05
params.method                   = "ls"
params.hsize                    = 50
params.celllines                = "K562,Nhlf,Hsmm,Gm12878,Hmec,Hepg2,Huvec,Nhek,H1hesc"

params.win                      = 10000
params.max                      = 0
params.dpms			= ""

params.type			= "differential"
params.pcutoff                  = 1e-5
params.niter                    = 25

params.dmrs                     = ""


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

    if (params.tool =="methylkey"){

       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --pdata pdata.txt --idat idat_rep/ --samples samples_names --groups samples_groups"
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--pdata                FILE                 Sample sheet"
       log.info "--idat                 FOLDER/FILE          Folder containing idat files (or file with beta values if platform=matrix)"
       log.info "--samples              STRING/INTEGER       Column name or column index for samples names"
       log.info "--groups               FILE                 Column name or column index for samples groups (eg : Case/Control)"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--html                 STRING               Name for output report (default=report.html)"
       log.info "--out                  STRING               Name for output directory (default=report_files)"
       log.info "--platform             STRING               IlluminaHumanMethylation450k, IlluminaHumanMethylationEPIC or matrix"
       log.info "--genome               STRING               Reference genome (default hg19)"
       log.info "--barcode              STRING/INTEGER       Column name or column index for array barcodes (default=barcode)"
       log.info "--pipeline             STRING               minfi or methylumi (default=minfi)"
       log.info "--normalize            STRING               Default is funnorm with minfi, bmiq with methylumi, use none to desactivate"
       log.info "--filters              FILE                 File containing list of probes to filter (one probe id by line)"
       log.info "--nalimit              DOUBLE               Maximum fraction of missing values for a probe"
       log.info "--missing              STRING               keep, mean (calculate mean of betas by group) or impute (with pamr)"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Flags:"
       log.info "--violin                                    Draw the violin plots of beta values by groups"
       log.info "--help                                      Display this message"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Cell Mixture Composition:"
       log.info "--cell houseman                             Houseman, only with minfi package, only for blood cells"
       log.info "--nbc                  INTEGER              Number of cell type (default=6)"
       log.info ""
       log.info ""
       exit 1
   }

   if (params.tool =="batchcorrection"){

       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --tool batchcorrection --meth methylkey/meth.rdata"
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--meth                 FILE                 meth.rdata file generated by methylkey"
       log.info "--variables            FILE                 variables to protect separeted by comma"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--correction           STRING               sva(default) or combat"
       log.info "--batch                STRING               With combat the variable name to correct"
       log.info "--what                 STRING               betas (correct original values, default) or mval (performe a new correction on the last generated mvalues)"
       log.info ""
       log.info ""
       exit 1

   }

   if (params.tool =="dmps"){

       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --tool dmps --meth methylkey/meth.rdata "
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--meth                 FILE                 meth.rdata file generated by methylkey"
       log.info "--variables            STRING               variables of your model"
       log.info "--case                 STRING               case"
       log.info "--control              STRING               control"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--fdr                  STRING               fdr cutoff for dmps (default=0.05)"
       log.info "--method               STRING               ls (default) or robust"
       log.info "--hsize                STRING               size of heatmap (default 50)"
       log.info "--model                STRING               custom model"
       log.info ""
       exit 1

   }

   if (params.tool =="comet"){
       
       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --tool comet --meth methylkey/meth.rdata --max 10"
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--meth                 FILE                 meth.rdata file generated by methylkey"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--win                  STRING               windows size (default=10000)"
       log.info "--max                  STRING               number of dmps to plot (default 0)"
       log.info "--dmps                 STRING               name of the dmps to plot"
       log.info ""
       exit 1
   }

   if (params.tool =="dmrcate"){
       
       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --tool dmrcate --meth methylkey/meth.rdata "
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--meth                 FILE                 meth.rdata file generated by methylkey"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--type                 STRING               differential (default) or variability"
       log.info "--fdr                  STRING               fdr cutoff for dmps (default=0.05)"
       log.info "--pcutoff              STRING               pvalue cutoff for dmrs (default=1e-5)"
       log.info ""
       exit 1

   }

   if (params.tool =="plotdmr"){
       
       log.info "--------------------------------------------------------"
       log.info "  USAGE                                                 "
       log.info "--------------------------------------------------------"
       log.info ""
       log.info "nextflow run iarcbioinfo/methylkey.nf --tool plotdmr --meth methylkey/meth.rdata --max 10"
       log.info ""
       log.info "Mandatory arguments:"
       log.info "--meth                 FILE                 meth.rdata file generated by methylkey"
       log.info ""
       log.info "--------------------------------------------------------"
       log.info "Optional arguments:"
       log.info "--max                  STRING               number of dmrs to plot (default 0)"
       log.info "--dmrs                 STRING               name of the dmrs to plot"
       log.info ""
       exit 1
   }


} else {

      /* Software information */
      if (params.tool =="methylkey"){
            log.info ""
            log.info "pdata                 = ${params.pdata}"
            log.info "idat                  = ${params.idat}"
            log.info "samples               = ${params.samples}"
            log.info "groups                = ${params.groups}"
            log.info "barcode               = ${params.barcode}"
            log.info "platform              = ${params.platform}"
            log.info "genome                = ${params.genome}"
            log.info "report files          = ${params.out}"
            log.info "pipeline              = ${params.pipeline}"
            log.info "normalize             = ${params.normalize}"
            log.info "filters               = ${params.filters}"
            log.info "NA cutoff             = ${params.nalimit}"
            log.info "missing values        = ${params.missing}"
            log.info "violin plot           = ${params.violin}"
            log.info ""
            log.info "CMC                   = ${params.cell}"
            log.info "number of cell type   = ${params.nbc}"
            log.info ""
      }

      if (params.tool =="batchcorrection"){
            log.info ""
            log.info "meth                  = ${params.meth}"
            log.info "variables             = ${params.variables}"
            log.info "correction            = ${params.correction}"
            log.info "batch                 = ${params.batch}"
            log.info "what                  = ${params.what}"
      }


      if (params.tool =="dmps"){
            log.info ""
            log.info "meth                  = ${params.meth}"
            log.info "variables             = ${params.variables}"
            log.info "model                 = ${params.model}"
            log.info "case                  = ${params.case}"
            log.info "control               = ${params.control}"
            log.info "fdr                   = ${params.fdr}"
            log.info "method                = ${params.method}"
            log.info "hsize                 = ${params.hsize}"
      }

      if (params.tool =="comet"){
            log.info ""
            log.info "meth                  = ${params.meth}"
            log.info "win                   = ${params.win}"
            log.info "max                   = ${params.max}"
            log.info "dmps                  = ${params.dmps}"
      }

      if (params.tool =="dmrcate"){
            log.info ""
            log.info "meth                  = ${params.meth}"
            log.info "type                  = ${params.type}"
            log.info "fdr                   = ${params.fdr}"
            log.info "pcutoff               = ${params.pcutoff}"
      }

      if (params.tool =="plotdmr"){
            log.info ""
            log.info "meth                  = ${params.meth}"
            log.info "max                   = ${params.max}"
            log.info "dmrs                  = ${params.dmrs}"
      }


}

if (params.tool == "methylkey"){

      pdata = file(params.pdata)
      idatDir = file(params.idat)

      flags=""
      if (params.violin){ flags = flags + " --violin" }

      process run_methylkey {

            publishDir params.out, mode: 'copy'

            input:

            file pdata
            file idatDir

            output:
            file("*") into report

            shell:
            '''
            methylkey.r --pdata !{pdata} --idat !{idatDir} --out . --samples !{params.samples} --groups !{params.groups} --pipeline !{params.pipeline} --nalimit !{params.nalimit} !{flags}
            '''
      }

}

if (params.tool == "batchcorrection"){

      meth = file(params.meth)

      process run_batchcorrection {

            publishDir params.out, mode: 'copy'

            input:
            file meth

            output:
            file("*") into report

            shell:
            '''
            batchcorrection.r --meth !{meth} --correction !{params.correction} --variables !{params.variables} --batch !{params.batch} --what !{params.what}
            '''
      }

}


if (params.tool == "dmps"){

      meth = file(params.meth)

      process run_dmps {

            publishDir params.out, mode: 'copy'

            input:
            file meth

            output:
            file("*") into report

            shell:
            '''
            dmps.r --meth !{meth} --model !{params.model} --variables !{params.variables} --case !{params.case} --control !{params.control} --fdr !{params.fdr} --method !{params.method} --hsize !{params.hsize}
            '''
      }

}


if (params.tool =="comet"){

      meth = file(params.meth)

      process run_comet {

            publishDir params.out, mode: 'copy'

            input:
            file meth

            output:
            file("*") into report

            shell:
            '''
            comet_plot.r --meth !{meth} --win !{params.win} --max !{params.max} --dmps !{params.dmps}
            '''
      }
}


if (params.tool =="dmrcate"){

      meth = file(params.meth)

      process run_dmrcate {

            publishDir params.out, mode: 'copy'

            input:
            file meth

            output:
            file("*") into report

            shell:
            '''
            dmrcate.r --meth !{meth} --type !{params.type} --fdr !{params.fdr} --pcutoff !{params.pcutoff}
            '''
      }
}

if (params.tool =="plotdmr"){

      meth = file(params.meth)
      dmrs = !{params.dmrs}
      if (dmrs != ""){ dmrs="--dmrs !{params.dmrs}"  } 

      process run_plotdmr {

            publishDir params.out, mode: 'copy'

            input:
            file meth

            output:
            file("*") into report

            shell:
            '''
            dmrs_plot.r --meth !{meth} --max !{params.max} !{dmrs}
            '''
      }
}









