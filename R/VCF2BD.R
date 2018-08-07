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
    x <- VCF53toBD_C(VCFFilename, outputFilename)
    tempDir <- strsplit(outputFilename, '/')[[1]]
    pattern <- paste(tempDir[length(tempDir)], ".tmp.", sep = "")
    if (length(tempDir) == 1)
      tempDir <- '.'
    else
      tempDir <- paste0(tempDir[1:length(tempDir) - 1], collapse = '/')
    currentDir <- getwd()
    setwd(tempDir)
    filesToDelete <- dir(pattern = pattern)
    file.remove(filesToDelete)
    setwd(currentDir)
    return (0)
  }
  return (VCF2BD_C(VCFFilename, infoFilename, outputFilename))
}
