BinaryDosage: Creates, Merges, and Reads Binary Dosage Files
================

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/USCbiostats/BinaryDosage?branch=UsingStandards&svg=true)](https://ci.appveyor.com/project/USCbiostats/BinaryDosage)
[![Travis build
status](https://travis-ci.org/USCbiostats/BinaryDosage.svg?branch=UsingStandards)](https://travis-ci.org/USCbiostats/BinaryDosage)
[![Codecov test
coverage](https://codecov.io/gh/USCbiostats/BinaryDosage/branch/UsingStandards/graph/badge.svg)](https://codecov.io/gh/USCbiostats/BinaryDosage?branch=UsingStandards)
<!-- badges: end -->

# Binary Dosage Files

### Introduction

Genotype imputation is an essential tool in genomics, enabling
association testing with markers not directly genotyped, increasing
statistical power, and facilitating data pooling between studies that
employ different genotyping platforms. Two commonly used software
packages for imputation are
[minimac](https://genome.sph.umich.edu/wiki/Minimac) and
[Impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html).
Furthermore, services such as the [Michigan Imputation
Server](https://imputationserver.sph.umich.edu/index.html) have made
genotype imputation much more accessible and streamlined.

While a number of software options are available for analyses of imputed
data (e.g. PLINK, EPACTS), fewer are available for Genomewide Gene x
Environment Interaction Scan (GWIS). Furthermore, data management tasks
such as parsing, subsetting, and merging, while manageable in smaller
studies, quickly become unyieldy and prohibitively slow with very large
samples sizes. We aim to address these limitations by converting
imputation outputs into a binary data format, BDose (Binary Dosage). The
benefits of a binary format are two fold - decreased hard drive storage
requirements (compared to a VCF file), and speed of parsing/analyses.
The BinaryDosage package contains functions to convert VCF and Impute2
formatted files into BDose format, along with functions to merge
samples.

For GWAS/GWIS analysis of BinaryDosage files, please refer to the
[**GxEScanR**](https://github.com/USCbiostats/GxEScanR) package.

### Description

##### Binary dosage data sets contain the following information:

  - Sample information
      - Family ID
      - Subject ID
  - SNP information
      - Chromosome number  
      - SNP ID  
      - Location in base pairs  
      - Reference allele  
      - Alternate allele  
  - Genetic information
      - Dosage values
      - Genotype probabilities, Pr(*g=0*), Pr(*g=1*), Pr(*g=2*)

There are 4 formats for a binary dosage data set. Data sets in formats
1, 2, and 3 have 3 files, a sample information file, a SNP information
file, and a genetic information file. Data sets in format 4 have just 1
file. This file contains all the information listed above and may
contain the following information.

  - Additional SNP information
      - Alternate allele frequency
      - Minor allele frequency
      - Average call rate
      - Imputation r squared
  - Merging information
      - Number of data sets merged
      - Sample size of each data set merged

### Functions

  - **vcftobd** - Converts a VCF file to abinary dosage data set
  - **gentobd** - Converts a GEN (impute2) file to a binary dosage data
    set
  - **bdmerge** - Merges multiple binary dosage data sets into a single
    data set
  - **getbdinfo** - Creates an R List containing information about a
    binary dosage data set (required for **getsnp** and **bdapply**)
  - **getvcfinfo** - Creates an R List containing information about a
    VCF file (required for **vcfapply**)
  - **getgeninfo** - Creates an R List containing information about a
    GEN file (required for **genapply**)
  - **bdapply** - Applies a function to the data for each SNP in a
    binary dosage file (requires list returned by **getbdinfo**)
  - **vcfapply** - Applies a function to the data for each SNP in a VCF
    file (requires list returned by **getvcfinfo**)
  - **genapply** - Applies a function to the data for each SNP in a GEN
    file (requires list returned by **getgeninfo**)
  - **getsnp** - Obtain genotype Dosages/Genotype Probabilities from a
    binary dosage file, outputs results to an R list

# Installation

1.  Install the [devtools](https://github.com/hadley/devtools) package
2.  Install the
    [BinaryDosage](https://github.com/USCbiostats/BinaryDosage) package
    directly from the USCbiostats repository on GitHub:

<!-- end list -->

``` r
remove.packages("BinaryDosage")
devtools::install_github("https://github.com/USCbiostats/BinaryDosage")

library(BinaryDosage)
```

# Usage

#### General Workflow

The general workflow for using binary dosage data sets is as follows:

  - Convert VCF or GEN files to a binary dosage data sets
      - Note: When converting a VCF file to a binary dosage data set,
        the information file associated with the vcf can be used to add
        addition data to the binary dosage data set
      - Note: When converting a GEN file to a binary dosage data set,
        the subject IDs can either be on the first line of the GEN file
        or in a separate sample file
  - Merge binary dosage datasets into a single data set
  - Apply a function to each SNP in the data set using bdapply
  - Extract SNPs for further analysis

#### Examples

##### Using VCF files

Example datasets *set1a.vcf* and *set1b.vcf* are representative of VCF
output files obtained from the Michigan Imputation Server. An
infromation file is also included. These files are in the extdata
directory of the BinaryDosage package. These files contain the
following:

| Study | N  | Number of SNPs |
| ----- | -- | -------------- |
| 1a    | 60 | 10             |
| 1b    | 40 | 10             |

###### Convert VCF files to the Binary Dosage Format and Merge into one data set

``` r
library(BinaryDosage)

# Get the name of the first VCF file
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
# Get the name of the first VCF information file
vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")

# Name of binary dosage file. Only need one because the default format is 4
bdfiles1a <- tempfile()
# Convert the file using the default format value of 4
vcftobd(vcffiles = c(vcf1afile, vcf1ainfo), bdfiles = bdfiles1a)

# Get the name of the second VCF file
vcf1bfile <- system.file("extdata", "set1b.vcf", package = "BinaryDosage")
# Get the name of the second VCF information file
vcf1binfo <- system.file("extdata", "set1b.info", package = "BinaryDosage")

# Name of binary dosage file. Only need one because the default format is 4
bdfiles1b <- tempfile()
# Convert the file using the default format value of 4
vcftobd(vcffiles = c(vcf1bfile, vcf1binfo), bdfiles = bdfiles1b)

# Merge data sets into one data set
mergebd1 <- tempfile()
bdmerge(mergefiles = mergebd1, bdfiles = c(bdfiles1a, bdfiles1b))
```

###### Apply function to all SNPs in data set

``` r
# Function to apply all SNPs in the merged data set This function calculates
# the alternate allele frequency
calculateaaf <- function(dosage, p0, p1, p2) {
    return(mean(dosage, na.rm = TRUE)/2)
}

# Get information about merged data set
mergebd1info <- getbdinfo(mergebd1)
aaf <- unlist(bdapply(mergebd1info, calculateaaf))
aafdf <- cbind(mergebd1info$snps, aaf)
aafdf
#>    chromosome location       snpid reference alternate      aaf
#> 1           1    10000 1:10000:C:A         C         A 0.352680
#> 2           1    11000 1:11000:T:C         T         C 0.013460
#> 3           1    12000 1:12000:T:C         T         C 0.240040
#> 4           1    13000 1:13000:T:C         T         C 0.337465
#> 5           1    14000 1:14000:G:C         G         C 0.190060
#> 6           1    15000 1:15000:A:C         A         C 0.562685
#> 7           1    16000 1:16000:G:A         G         A 0.456940
#> 8           1    17000 1:17000:C:A         C         A 0.457760
#> 9           1    18000 1:18000:C:G         C         G 0.259105
#> 10          1    19000 1:19000:T:G         T         G 0.243095

# Get the dosage values for the 6th SNP
snp <- getsnp(mergebd1info, 6)
snpdf <- data.frame(dosage = unlist(snp), row.names = mergebd1info$samples$sid)
head(snpdf)
#>    dosage
#> I1  1.000
#> I2  1.849
#> I3  1.000
#> I4  2.000
#> I5  1.046
#> I6  1.915
```
