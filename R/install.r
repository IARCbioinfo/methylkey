if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
install.packages("R.utils")

BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("wateRmelon")

install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")