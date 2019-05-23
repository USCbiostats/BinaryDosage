#' Function to read SNP values
#'
#' Function to read SNP values from VCF or Binary Dosage files.
#'
#' @param fileInfo
#' List returned from GetVCFInfo or GetBDoseInfo with value
#' index set to TRUE
#' @param SNPs
#' List of SNPs to read from the file. Can be either a character
#' or integer vector. Routine is much faster if a sorted integer
#' vector is provided.
#' @param subjects
#' List of subjects to read SNP values for. This can be either
#' a character or integer vector. If no value is provided SNP
#' values are read for all subjects
#' @param geneProb
#' Indicator if genetic probabilities are to be returned. The
#' default value is TRUE
#'
#' @return
#' A matrix with the dosage and genetic probabilities in the
#' columns and each row corresponds to one subject.
#' @export
#'
#' @examples
#' # Under construction
GetSNPValues <- function(fileInfo, SNPs, subjects, geneProb = TRUE) {
  if (missing(fileInfo) == TRUE)
    stop("No genetic file information specified")
  if ("genetic-file-info" %in% class(fileInfo) == FALSE)
    stop("fileInfo is not of class \"genetic-file-info\"")
  if (missing(SNPs) == TRUE)
    stop("No SNPs specified")
#  if (is.character(SNPs)) {
#
#  }
  return (0)
}
