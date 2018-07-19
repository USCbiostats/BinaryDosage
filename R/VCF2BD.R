#' Function to creae a binary dosage file from a VCF file
#'
#' Function to creae a binary dosage file from a VCF file
#'
#' @param VCFFilename
#' VCF file name
#' @param outputFilename
#' Name of output file
#' @param infoFilename
#' Name of info file associated with VCF file. Can be blank.
#' Default value = ""
#' @param calculateMAF
#' Indicator on whether to calculate MAF based on the sample.
#' Default value = TRUE
#' @return
#' 0 - success
#' 1 - failure
#' @export
VCF2BD <- function(VCFFilename, outputFilename, infoFilename = "", calculateMAF = TRUE) {
  if (infoFilename == "") {
    return (VCF53toBD_C(VCFFilename, outputFilename))
  }
  return (VCF2BD_C(VCFFilename, infoFilename, outputFilename))
}
