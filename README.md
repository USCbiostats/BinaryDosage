BinaryDosage: Convert Imputed SNP data to a Binary Dosage Format
================

BDose Format
============

### Introduction

Genotype imputation is an essential tool in genomics, enabling association testing with markers not directly genotyped, increasing statistical power, and facilitating data pooling between studies that employ different genotyping platforms. Two commonly used software packages for imputation are [minimac](https://genome.sph.umich.edu/wiki/Minimac) and [Impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html). Furthermore, services such as the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html) have made genotype imputation much more accessible and streamlined.

While a number of software options are available for analyses of imputed data (e.g. PLINK, EPACTS), fewer are available for Genomewide Gene x Environment Interaction Scan (GWIS). Furthermore, data management tasks such as parsing, subsetting, and merging, while manageable in smaller studies, quickly become unyieldy and prohibitively slow with very large samples sizes. We aim to address these limitations by converting imputation outputs into a binary data format, BDose (Binary Dosage). The benefits of a binary format are two fold - decreased hard drive storage requirements (compared to a VCF file), and speed of parsing/analyses. The BinaryDosage package contains functions to convert VCF and Impute2 formatted files into BDose format, along with functions to merge and subset samples/SNPs.

For GWAS/GWIS analysis of BDose files, please refer to the [**GxEScanR**](https://github.com/USCbiostats/GxEScanR) package.

### Description

##### BDose files contain at minimum the following information:

-   Subject IDs
-   SNP information
-   Chromosome number
-   Location in base pairs
-   Reference Allele
-   Alternate Allele
-   Dosage values

##### BDose files may contain any of the following information:

-   Family IDs
-   SNP information
-   SNP name (rs number)
-   Minor allele indicator (reference of alternate allele)
-   Minor allele frequency
-   Average call rate
-   r2
-   Genetic probabilities, Pr(g=0), Pr(g=1), Pr(g=2)

##### MERGED BDose files contain the following information:

-   Number of files combined
-   Number of subjects in each VCF file

### Functions

-   **VCFtoBD** - Convert VCF/Impute2 files to BDose format
-   **MergeBD** - Merge BDose files
-   (Currently retains only intersection of markers based on *Chromosome:BasePair\_RefAllele\_AltAllele*)
-   **GetBDoseInfo** - Create R List containing BDose file attributes (required for **GetSNPValues**)
-   **GetSNPValues** - Obtain genotype Dosages/Genotype Probabilities, outputs results to an R Data Frame

Installation
============

1.  Install the [devtools](https://github.com/hadley/devtools) package
2.  Install the [BinaryDosage](https://github.com/USCbiostats/BinaryDosage) package directly from the USCbiostats repository on GitHub:

``` r
library(devtools)
install_github("USCbiostats/BinaryDosage")

library(BinaryDosage)
```

Usage
=====

#### General Workflow

##### [**BinaryDosage**](https://github.com/USCbiostats/BinaryDosage)

-   Convert VCF to BDose

##### [**GxEScanR**](https://github.com/USCbiostats/GxEScanR)

-   Use *ReadBDInfo()* to create an R list of BDose file attributes
-   Note: This is the only way to see what is in a BDose file
-   *read.table* phenotype data (*R*)
-   Note: order of covariates matter, refer to package documentation
-   Run GxEScanR to conduct GWAS/GWIS

#### Example

Example datasets *batch1.vcf.gz* and *batch2.vcf.gz* are representative of VCF output files obtained from the Michigan Imputation Server, and contain the following:

| Study | N   | Number of SNPs |
|-------|-----|----------------|
| B1    | 100 | 500            |
| B2    | 100 | 500            |

###### Convert VCF to BDose Format

``` r
library(BinaryDosage)

# Currently, must extract compressed files (vcf.gz support TBD)
# system('gzcat batch1.vcf.gz > batch1.vcf') system('gzcat batch2.vcf.gz >
# batch2.vcf')

# Convert VCF Batch 1 and Batch 2 files:
VCFtoBD("batch1.vcf", "batch1.bdose")
VCFtoBD("batch2.vcf", "batch2.bdose")
```

<!-- ``` {r, eval = T, message = F, collapse = T} -->
<!-- bdose1 <- GetBinaryDosageInfo("./example/batch1.bdose") -->
<!-- -->
<!-- bdose1$NumSamples -->
<!-- bdose1$NumSNPs -->
<!-- ``` -->
###### Merge BDose files

``` r
library(BinaryDosage)
filesToMerge <- c("batch1.bdose", "batch2.bdose")
MergeBD("merged.bdose", filesToMerge)
```

###### Verify contents

``` r
pooled <- GetBDoseInfo("merged.bdose")
pooled$NumSamples
#> [1] 200
pooled$NumSNPs
#> [1] 500
```

###### GWAS/GWIS - See [**GxEScanR**](https://github.com/USCbiostats/GxEScanR)
