#' Function to read in SNPs from a binary dosage file
#'
#' Function to read in SNPs from a binary dosage file
#'
#' @param bdInfo
#' Information about the binary dosage file return from
#' GetBDoseInfo
#' @param SNPs
#' Vector of either SNP names or indices in binary dosage
#' file to extract values for
#' @param Subjects
#' Vector of either Subject IDs of indices in binary dosage
#' file to extract values for
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetSNPValues <- function(bdInfo, SNPs, Subjects) {
  if (missing(bdInfo))
    return (NULL)
  if (missing(SNPs)) {
    return (NULL)
  } else {
    if (is.vector(SNPs, "character")) {
      snpVec <- match(SNPs, bdInfo$SNPs$SNPID)
    } else if (is.vector(SNPs, "integer")) {
      snpVec <- SNPs
    } else {
      return (NULL)
    }
  }
  if (missing(Subjects)) {
    subVec <- 1L:nrow(bdInfo$Samples)
  } else {
    if (is.vector(Subjects, "character")) {
      subVec <- match(Subjects, bdInfo$Samples$SID)
    } else if (is.vector(Subjects, "integer")) {
        subVec <- Subjects
    } else {
        return (NULL)
    }
  }
  if (anyNA(snpVec))
    return (NULL)
  if (anyNA(subVec))
    return (NULL)

  if(bdInfo$version == 2)
    valueMatrix <- matrix(0, length(subVec), 4 * length(snpVec))
  else
    valueMatrix <- matrix(0, length(subVec), length(snpVec))

  GetSNPValuesC(bdInfo$filename, subVec, snpVec, bdInfo$Indices, valueMatrix)

  return (list(snpVec, subVec, valueMatrix))
}
