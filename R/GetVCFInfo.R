#' Function to read in information about a vcf file
#'
#' Function to read in information about a vcf file.
#' This information is used when reading the dosage and genetic
#' probabilities from the file
#'
#' @param vcfFile
#' Name of the vcf file
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetVCFInfo <- function(vcfFile, reserve) {
  if (missing(vcfFile))
    return (NULL)
  if (missing(reserve))
    reserve <- 0
  p1 <- BinaryDosage:::GetVCFHeaderC(vcfFile)
  p2 <- BinaryDosage:::GetVCFSNPInfoC(p1$filename, p1$StartData, reserve)
  return (c(p1,p2))
}
