# methylkey
## pipeline for bisulfite data analysis

## Description
methylkey is a pipeline for 450k and 850k array analysis using Minfi; Methylumi, Comet, Bumphunter and DMRcate packages.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. R packages: 

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

# Loading data and quality control

## Input
  | Type      | Description           |
  |-----------|-----------------------|
  | --pdata   | Sample sheet          |
  | --idat    | Idat files directory  |

## Parameters

  * #### Mandatory
| Name        | Example value        | Description                                                     |
|-------------|----------------------|-----------------------------------------------------------------|
| --samples   |  Sample_name         | Name of the column containing sample's name in the Sample sheet |
| --groups    |  Sample_Group,Gender | Name of the columns containing groups in the Sample sheet       |
| --barcode   |  barcode             | Name of the column containing barcodes in the Sample sheet      | 

  * #### Optional
| Name        | Default value            | Description                                |
|-------------|--------------------------|--------------------------------------------|
| --pipeline  | minfi                    | minfi (default) or methylumi               |
| --filters   | Crossreactive_probes.csv | file(s) containing probes id to remove     |
| --nalimit   | 0.2                      | maximum fraction of Na values for a probe  |
| --missing   | mean                     | mean, input or keep                        |

  * #### Flags

Flags are special parameters without value.

| Name      | Description                       |
|-----------|-----------------------------------|
| --help    | Display help                      |
| --violin  | report violin plot of beta values |

## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --pdata pdata --idat idat --samples Sample_Name --groups Sample_Group,Gender --filters ../data/Crossreactive_probes.csv
  ```
  
## Output
  | Type      | Description                           |
  |-----------|---------------------------------------|
  | --report  | report folder name (default methylkey)|
  
  
# Batch Correction

## Input
  | Type        | Description                          |
  |-------------|--------------------------------------|
  | --meth      | meth.rdata file in report folder     |

## Parameters

  * #### Mandatory
| Name        | Example value        | Description                                 |
|-------------|----------------------|---------------------------------------------|
| --variables |  Sample_Group,Gender | variables to protect ( separeted by comma ) |

  * #### Optional
| Name         | Default value  | Description                               |
|--------------|----------------|-------------------------------------------|
| --correction | sva            | sva(default) or combat                    |
| --batch      |                | With combat the variable name to correct  |

## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --tool batchcorrection --meth methylkey/meth.rdata --variables Sample_Group,Gender
  ```


# DMPs

## Input
  | Type        | Description                          |
  |-------------|--------------------------------------|
  | --meth      | meth.rdata file in report folder     |

## Parameters

  * #### Mandatory
| Name        | Example value        | Description                                 |
|-------------|----------------------|---------------------------------------------|
| --variables |  Sample_Group,Gender | variables to protect ( separeted by comma ) |
| --case      |  Case                | case group                                  |
| --control   |  Control             | control group                               |

  * #### Optional
| Name         | Default value  | Description                               |
|--------------|----------------|-------------------------------------------|
| --fdr        | 0.05           | fdr cutoof for dmps                       |
| --method     | ls             | ls or robust                              |
| --hsize      | 50             | size of the heatmap                       |


## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --tool dmps --meth methylkey/meth.rdata --variables Sample_Group,Gender
 --case case --control control  
 ```

# CoMet

## Input
  | Type        | Description                          |
  |-------------|--------------------------------------|
  | --meth      | meth.rdata file in report folder     |

## Parameters

  * #### Mandatory
| Name        | Example value           | Description                                 |
|-------------|-------------------------|---------------------------------------------|
| --max       |  10                     | print the first 10 dmps                     |
| --dmps      |  cg06955954,cg23899408  | print dmps by ids                           |

use --max, --dmps or both.

  * #### Optional
| Name         | Default value  | Description                               |
|--------------|----------------|-------------------------------------------|
| --win        | 10000          | windows size in base pair                 |

## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --tool comet --meth methylkey/meth.rdata --max 10
  ```
  
# DMRcate

## Input
  | Type        | Description                          |
  |-------------|--------------------------------------|
  | --meth      | meth.rdata file in report folder     |

## Parameters

  * #### Optional
| Name         | Default value  | Description                                 |
|--------------|----------------|---------------------------------------------|
| --type       | differential   | analysis type(differential or variability ) |
| --fdr        | 0.05           | fdr cutoff for dmps                         |
| --pcutoff    | 1e-5           | pvalue cutoff for dmrs                      |

## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --tool dmrcate --meth methylkey/meth.rdata
  ```
  
# plot DMRs

## Input
  | Type        | Description                          |
  |-------------|--------------------------------------|
  | --meth      | meth.rdata file in report folder     |

## Parameters

 * #### Mandatory
| Name        | Example value           | Description                                 |
|-------------|-------------------------|---------------------------------------------|
| --max       |  10                     | print the first 10 dmrs                     |

use --max, --dmps or both.

## Usage
  ```
  nextflow run iarcbioinfo/methylkey.nf --tool plotdmr --meth methylkey/meth.rdata --max 10
  ```


