# methylkey: DNA Methylation Analysis from Illumina Arrays

[![DOI](https://zenodo.org/badge/139563728.svg)](https://zenodo.org/badge/latestdoi/139563728)

`methylkey` is a Bioconductor package for comprehensive analysis of DNA methylation data from Illumina methylation arrays. It supports multiple platforms (27k, 450k, EPIC, EPIC+, EPICv2, and mouse arrays) with preprocessing pipelines using both `sesame` and `minfi` packages.

## Features

- **Flexible Preprocessing**: Choose between sesame or minfi pipelines
- **Quality Control**: Comprehensive QC plots and statistics
- **Differential Methylation**: DMP and DMR detection using multiple methods
- **Platform Detection**: Automatic detection of array platform
- **Sex Inference**: Automatic inference of sample sex
- **Epigenetic Age**: Multiple clock model support
- **Cell Type Estimation**: Blood immune cell composition estimation
- **Standardized Output**: SummarizedExperiment-based objects for downstream analysis

## Installation

### From GitHub

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_github("IARCbioinfo/methylkey")
library(methylkey)
```

### Requirements

- R >= 4.2.1
- Bioconductor packages (automatically installed)
- For minfi pipeline: `minfi`, `wateRmelon`
- For sesame pipeline: `sesame`, `sesameData`

## Quick Start

### 1. Prepare Sample Sheet

Create a tab-separated file with sample information:

```
Sample_Name  Basename                  Group   Age
Sample_001   203021070069_R03C01       Ctrl    45
Sample_002   203021070069_R04C01       Case    50
Sample_003   203021070069_R05C01       Ctrl    48
```

**Required columns:**
- `Basename` or `barcode`: Sentrix ID and position (e.g., `203021070069_R03C01`)
- At least one additional column for grouping or clinical variables

Load the sample sheet in R:

```r
sampleSheet <- readr::read_delim("path/to/sampleSheet.txt")
```

### 2. Preprocessing with sesame

For mouse arrays or when sesame is preferred:

```r
library(methylkey)

# Load IDAT files and preprocess with sesame
meth <- sesame2Betas(
  idat = "path/to/idat/directory",
  sampleSheet = sampleSheet,
  prep = "QCDPB",
  na = 0.2,
  ncore = 4
)

# View object information
meth
colData(meth)
get_betas(meth)[1:5, 1:5]
```

### 3. Preprocessing with minfi

For 450k and EPIC arrays:

```r
meth <- minfi2Betas(
  idat = "path/to/idat/directory",
  sampleSheet = sampleSheet,
  na = 0.2,
  compositeCellType = "Blood"  # Optional: estimate cell counts
)
```

### 4. Quality Control

```r
# Access QC statistics
qcs <- metadata(meth)$qcs
head(qcs)

# Generate QC plots
plot_Detection(meth)
plot_NA(meth)
```

### 5. Differential Methylation Analysis

```r
# Convert to M-values and perform analysis
mvals <- get_mvals(meth, grp = "group")

# Run differential methylation analysis
mrs <- methyldiff(
  se = meth,
  model = "~group",
  intercept = "Ctrl",
  method = "ls",
  fdr = 0.05,
  pcutoff = 0.2
)

# Extract results
dmps <- get_dmps(mrs)
dmrs <- get_dmrs(mrs)
```

## Supported Array Platforms

| Platform | Probes | Support |
|----------|--------|---------|
| IlluminaHumanMethylation27k | 27,578 | ✓ |
| IlluminaHumanMethylation450k | 486,427 | ✓ |
| IlluminaHumanMethylationEPIC | 866,553 | ✓ |
| IlluminaHumanMethylationEPIC+ | 865,546 | ✓ |
| IlluminaHumanMethylationEPICv2 | 937,690 | ✓ |
| IlluminaMouseMethylation285k | 296,070 | ✓ |

## Data Objects

### Betas Object

Container for beta values (0-1 scale):
```r
# Access
betas <- get_betas(meth)
samples <- colData(meth)
platform <- metadata(meth)$plateform
```

### Mvals Object

Container for M-values (log2 scale) for statistical analysis:
```r
mvals <- get_mvals(betas)
mvals_matrix <- assays(mvals)$mvals
```

## Documentation

- Full function documentation: `?methylkey` or `help(package="methylkey")`
- Vignettes: `browseVignettes("methylkey")`
- Examples: See `examples/` directory for complete workflows

## Workflow Summary

```
IDAT Files
    ↓
[sesame2Betas or minfi2Betas]
    ↓
Betas Object
    ↓
[QC & Visualization]
    ↓
[get_mvals]
    ↓
Mvals Object
    ↓
[methyldiff]
    ↓
MethylResultSet (DMPs & DMRs)
```

## Common Issues

**"Package sesame not found"**
```r
BiocManager::install("sesame")
```

**"Basename not in sampleSheet"**
Ensure your sample sheet has a `Basename` column with barcode information.

**"Barcodes in sampleSheet not found in betas"**
Verify that barcodes in the sample sheet match IDAT file names exactly.

## Citation

If you use methylkey in your research, please cite:

```bibtex
@software{methylkey2026,
  author = {Cahais, Vincent},
  year = {2026},
  title = {methylkey: DNA methylation analysis from Illumina arrays},
  url = {https://github.com/IARCbioinfo/methylkey}
}
```

## References

- Aryee MJ, et al. (2014). Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays. Bioinformatics, 30(10), 1363–1369.
- Zhou W, et al. (2018). SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. NAR, 46(20), e119.
- Hansen KD, et al. (2016). Reference-based correction of interstrand dye bias in two-color microarrays. Genome Biology, 17(1), 37.

## License

GPL-3.0

## Contributing

Contributions are welcome! Please open an issue or pull request on GitHub.

## Contact

For questions and issues: [GitHub Issues](https://github.com/IARCbioinfo/methylkey/issues)
