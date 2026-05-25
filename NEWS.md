CHANGELOG - methylkey package improvements
Version 0.99.0 (2026-04-29)

## Major Changes

### Package Structure and Metadata
- Created comprehensive DESCRIPTION file with all dependencies properly listed
- Updated NAMESPACE with roxygen imports and exports
- Enhanced methylkey-package.r with complete documentation and import statements
- Updated .Rbuildignore to exclude development artifacts
- Added package website metadata (URL, BugReports)

### Bug Fixes and Code Quality

#### sesame2Betas():
- Removed `saveRDS(sdfs, "sdfs.rds")` debug artifact
- Replaced `require()` with `requireNamespace()` for proper error handling
- Added namespace qualification to all sesame function calls (sesame::)
- Improved error messages and logging
- Added handling for optional Clock_models parameter
- Added metadata for preprocessing method
- Improved message output for pipeline tracking

#### minfi2Betas():
- Removed all debug messages (message("001"), etc.)
- Replaced `require()` with `requireNamespace()`
- Improved platform-specific cell count estimation logic
- Added proper error handling for missing optional packages
- Enhanced immune score calculation with null checks
- Added platform detection to metadata
- Better handling of extended cell type models

#### newBetas():
- Fixed error message handling (was using cat() in assert_that)
- Improved validation with better error messages
- Added numeric validation for betas matrix
- Enhanced barcode validation with informative messages
- Better handling of missing samples in data

#### formatSampleSheet():
- Enhanced robustness with explicit assertions
- Added check for empty data frames
- Improved barcode/basename column handling
- Better handling of sentrix column separation
- Added proper null check for column operations
- Improved documentation

#### toGeoSubmission() for Betas:
- Completely rewritten to work correctly
- Uses proper getBetas() extraction
- Generates proper TSV output
- Includes sample names in headers
- Added validation and error handling

#### Other Improvements:
- Imported all sesame functions that were previously using bare calls
- Added BiocParallel to imports for parallel processing
- Fixed several implicit function dependencies

### Testing Infrastructure

#### New test suite (tests/testthat/):
- test_sampleSheet.R: Tests for formatSampleSheet()
  - Basic input handling
  - Column renaming
  - Error handling for invalid inputs
  
- test_platform.R: Tests for getPlateform()
  - EPIC detection
  - 450k detection
  - Mouse 285k detection
  - Unknown platform handling
  - NULL input handling

- test_newBetas.R: Tests for newBetas()
  - Valid Betas object creation
  - Input validation
  - Barcode requirement
  - NA parameter validation

- test_missingValues.R: Tests for CpGNAexcl()
  - NA probe identification
  - Nalimit parameter behavior
  - Empty NA handling

### Documentation

#### README.md:
- Complete rewrite with modern format
- Installation instructions
- Quick start guide
- Sample sheet preparation details
- Code examples for both sesame and minfi pipelines
- QC and differential analysis examples
- Supported platforms table
- Data object documentation
- Troubleshooting section
- References and citations

#### Package-level Documentation (methylkey-package.r):
- Added comprehensive roxygen documentation
- Listed all main features
- Provided workflow overview
- Clear import specifications

### Dependency Management

Updated and clarified dependencies:
- **Imports**: Added all required functions explicitly
  - sesame, sesameData, minfi, wateRmelon
  - Data processing: dplyr, tidyr, readr, purrr
  - Visualization: ggplot2, ggpubr, ggvenn, Gviz
  - Statistics: rstatix, randomForest, limma, sva
  - Core: SummarizedExperiment, GenomicRanges, S4Vectors
  - Utilities: assertthat, DescTools

- **Suggests**: Added optional packages
  - knitr, rmarkdown, testthat
  - FlowSorted.Blood.EPIC (for cell estimation)
  - BiocStyle, pkgdown
  - BSgenome packages (for genomic operations)

### Removed/Deprecated

- Removed hard-coded `require()` calls in favor of `requireNamespace()`
- Removed debug messages and artifacts
- Removed broken toGeoSubmission method for RGChannelSet (if not functional)

## Known Limitations

- Clock model support in sesame2Betas requires external model files
- Cell count estimation requires external reference packages
- Some mouse-specific functions may not work with human arrays

## Next Steps (Phase 6-7)

- Add vignettes for preprocessing workflows
- Implement methylkey_report() function with Rmd template
- Set up GitHub Actions for CI/CD
- Implement pkgdown website
- Add more comprehensive integration tests
- Optimize memory usage for large datasets

## Installation Instructions for Testing

```r
# From directory
devtools::load_all()
devtools::test()
devtools::check()
```

## Contact

For detailed implementation notes and technical discussions, see ROADMAP.txt
