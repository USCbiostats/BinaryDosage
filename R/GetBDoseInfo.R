#' Function to read in information about binary dosage file
#'
#' Function to read in information about binary dosage file.
#' This information is used when read the dosage and genetic
#' probabilities from the file
#'
#' @param bdFile
#' Name of the binary dosage file
#' @param famFile
#' Name of file with subject data. This file is in the format
#' used by plink .fam files. Not needed for format 4 or greater
#' binary dosage files.
#' @param mapFile
#' Name of file with SNP information. This file is in the format
#' used by plink for extended map file files. Not needed for
#' format 4 or greater binary dosage files.
#' @param index
#' Index the SNPs for faster reading.
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetBDoseInfo <- function(bdFile, famFile = "", mapFile = "", index = FALSE) {
  if (missing(bdFile))
    return (NULL)
  if (index == TRUE)
    x = 1L
  else
    x = 0L
  return (GetBinaryDosageInfoC(bdFile, famFile, mapFile, x))
}
