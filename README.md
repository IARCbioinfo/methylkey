# methylkey
## pipeline for bisulfite data analysis

## Description
methylkey is a pipeline for 450k and 850k array analysis using Minfi; Methylumi, Comet, Bumphunter and DMRcate packages.

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. R packages: 

You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.

## Input
  | Type      | Description           |
  |-----------|-----------------------|
  | --pdata   | Sample sheet          |
  | --idat    | Idat files directory  |

  Specify the test files location

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
  | Type      | Description        |
  |-----------|--------------------|
  | report    | report folder name |

