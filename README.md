# Methylkey
#### 05 Janvier 2021

1. Installation
2. Prepare SampleSheet
3. Preprocessing illumina array
  3.1 Minfi (450k and EPIC)
  3.2 Sesame (285k Mouse)

4. Plot QC report
5. Differential methylation analysis

## 1. Installation :

```bash
install.packages("devtools") 
library(devtools)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

devtools::install_github("IARCbioinfo/methylkey@v1.0")
library(methylkey)
```

## 2. Prepare sampleSheet

The sampleSheet is a tabular file with samples as lines and variables as columns. This file must include : 
  * samples names ( uniques values )
  * Barcode (eg: 203021070069_R03C01) : this code is the link between samples and idat files. The fist number (203021070069) is the sentrix id, while (R03C01) is the sentrix position. Each sample is associated with two idat files, one for the red channel (*203021070069_R03C01_red.idat*) and one for the green channel (*203021070069_R03C01_green.idat*).
  * Any variables that describe samples : Treatment, Status, Age, CellType, ...


## 3. Preprocessing illumina array

### 3.1 Minfi (450k and EPIC)

[The minfi User’s Guide](http://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html)

[Example](examples/dataloader.minfi.Rmd)

### 3.2 Sesame (285k Mouse)

[SeSAMe Basic Usage](https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html)

[Example](examples/dataloader.minfi.Rmd)




## Bibliography

Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA (2014). “Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays.” Bioinformatics, 30(10), 1363–1369. doi: 10.1093/bioinformatics/btu049.

Maksimovic J, Gordon L, Oshlack A (2012). “SWAN: Subset quantile Within-Array Normalization for Illumina Infinium HumanMethylation450 BeadChips.” Genome Biology, 13(6), R44. doi: 10.1186/gb-2012-13-6-r44.

Fortin J, Labbe A, Lemire M, Zanke BW, Hudson TJ, Fertig EJ, Greenwood CM, Hansen KD (2014). “Functional normalization of 450k methylation array data improves replication in large cancer studies.” Genome Biology, 15(12), 503. doi: 10.1186/s13059-014-0503-2.

Triche TJ, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD (2013). “Low-level processing of Illumina Infinium DNA Methylation BeadArrays.” Nucleic Acids Research, 41(7), e90. doi: 10.1093/nar/gkt090.

Fortin J, Hansen KD (2015). “Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data.” Genome Biology, 16, 180. doi: 10.1186/s13059-015-0741-y.

Andrews SV, Ladd-Acosta C, Feinberg AP, Hansen KD, Fallin MD (2016). “'Gap hunting' to characterize clustered probe signals in Illumina methylation array data.” Epigenetics & Chromatin, 9, 56. doi: 10.1186/s13072-016-0107-z.

Fortin J, Triche TJ, Hansen KD (2017). “Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi.” Bioinformatics, 33(4). doi: 10.1093/bioinformatics/btw691.

Zhou W, Triche TJ, Laird PW, Shen H (2018). “SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions.” Nucleic Acids Research, gky691. doi: 10.1093/nar/gky691.
