#' Function to read in information about a vcf file
#'
#' Function to read in information about a vcf file.
#' This information is used when reading the dosage and genetic
#' probabilities from the file
#'
#' @param vcfFile
#' Name of the vcf file
#' @param bdFile
#' Name of the binary dosage file
#' @param famFile
#' Name of the family file, only used with formats 1, 2, and 3
#' @param mapFile
#' Name of the map file, only used with formats 1, 2, and 3
#' @param format
#' Format of binary dosage file - default 4, can be 1, 2, or 3
#' @param version
#' Version of the format to use - default 2, can be 1
#' @param batchSize
#' Size of batches of SNPs to read in before processing
#' @return
#' 0 - Success
#' 1 - Failure
#' @export
VCFtoBD <- function(vcfFile, bdFile, famFile = "", mapFile = "", format = 4, version = 2, batchSize = 100) {
  if (missing(vcfFile))
    return (NULL)
  if (missing(bdFile))
    return (NULL)
  vcfi <- GetVCFInfo(vcfFile)
  if (length(vcfi) == 0)
    return (NULL)
  return (BDConvertVCFC(vcfi, bdFile, famFile, mapFile, format, version, batchSize))
}
