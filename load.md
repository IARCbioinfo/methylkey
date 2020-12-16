
# Loading data and quality control

This will load your data into methylkey from idat files.


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
