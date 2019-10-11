#!/usr/bin/env nextflow
echo true

params.path="~/git/methylkey/bin/"
params.bad="bad_samples.txt"
params.bed="None"
params.analysisPath="/mojo/EPIMARK/analyse_20191009/"

/* Build all analysis setups */
idats=Channel.from([ ["cpgislands", "$params.analysisPath/betas_cpgislands.csv.gz","None"],
                     ["geneBody","$params.analysisPath/betas_genes.csv.gz","None"],
                     ["promoter","$params.analysisPath/betas_promoters.csv.gz","None"] ])      
tissues=Channel.from([ ["tissue","$params.analysisPath/pdata.csv"] ])
    .spread(idats)

idats=Channel.from([ ["cpgislands", "$params.analysisPath/betas_BC_cpgislands.csv.gz","bad_samples.txt"],
                     ["geneBody","$params.analysisPath/betas_BC_genes.csv.gz","bad_samples.txt"],
                     ["promoter","$params.analysisPath/betas_BC_promoters.csv.gz","bad_samples.txt"] ])
blood=Channel.from([ ["blood","$params.analysisPath/pdata_BC.csv"] ])
    .spread(idats)

tissues.concat( blood ).set { step1 }

models=[ [ 'crude', '3' ], ['paired','c(3,4)'] ]

/*methylkey QC for each idat*/
process methylkey {

  publishDir 'report', mode: 'symlink'

  input: 
    set val(material), val(pdata), val(region), val(idat), val(bad) from step1

  output:
    set val(material), val(pdata), val(region), val(idat) into step2      
    file 'qc_report*' into qc_report
    file 'qc*html' into qc_html

  script:
    """
    cp $params.path/methylkey_rnbeads.Rmd qc_${material}_${region}.Rmd;
    Rscript -e "rmarkdown::render('qc_${material}_${region}.Rmd', params=list( path='$params.path' , \
                        pdata='$pdata' , idat='$idat' , genome='hg38' ,  samples=1 , groups=3 , nalimit=0.05 , \
                        missing='mean' , out='qc_report_${material}_${region}' , violin=FALSE , badsamples='$bad' ))"
    """
}


/*SVA for each idat using two models*/
process sva {

  publishDir 'report', mode: 'symlink'

  input:
    file report from qc_report
    set val(material), val(pdata), val(region), val(idat) from step2
    each correction from models

  output:
    set val(material), val(pdata), val(region), val(idat), val(correction) into step3
    file 'sva_report*' into sva_report
    file 'sva*html' into sva_html

  script:
    index="sva_" + material + "_" + region + "_" + correction[0] + ".Rmd"
    rep="sva_report_" + material + "_" + region + "_" + correction[0]
    corr=correction[1]
    """
    cp $params.path/methylkey_batchcorrection.Rmd $index
    Rscript -e "rmarkdown::render('$index', params=list(  path='$params.path' , \
                        meth='$report/meth.rdata' , variables=$corr , out='$rep' , \
                        correction='sva' , nalimit='0.05' , missing='mean'))"
    """
}

/*DMRS for each sva using two models*/
process dmrs {

  publishDir 'report', mode: 'symlink'

  input: 
    file report from sva_report
    set val(material), val(pdata), val(region), val(idat), val(correction) from step3
    each model from models

  output:
    file 'dmrs_report*' into dmrs_report
    file 'dmrs*html' into dmrs_html

  script:
    index="dmrs_" + material + "_" + region + "_" + correction[0] + "_" + model[0] + ".Rmd"
    rep="dmrs_report_" + material + "_" + region + "_" + correction[0] + "_" + model[0]
    mod=model[1]
    """
    echo "what:$material\nregion:$region\nsva:$correction\nmodel:$model\n" > log
    cp $params.path/methylkey_dmps.Rmd $index
    Rscript -e "rmarkdown::render('$index', params=list(  path='$params.path' , genome='hg38' , \
                        meth='$report/meth.rdata' , variables=$mod , case='TT' , control='NT' , \
                        method='ls' , fdr=0.05 , hsize=50 , logistic=TRUE , bed='$params.bed' , out='$rep'))"
    """
}