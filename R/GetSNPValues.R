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

  if (GetSNPValuesC(bdInfo$filename, subVec, snpVec, bdInfo$Indices, valueMatrix) == 1)
    return (NULL)

  if (bdInfo$version == 2)
    colnames(valueMatrix) <- paste0(rep(bdInfo$SNPs$SNPID[snpVec], each = 4), rep(c("_dose", "_p0", "_p1", "_p2"), length(snpVec)))
  else
    colnames(valueMatrix) <- bdInfo$SNPs$SNPID[snpVec]
  rownames(valueMatrix) <- bdInfo$Samples$SID[subVec]

  return (valueMatrix)
}
