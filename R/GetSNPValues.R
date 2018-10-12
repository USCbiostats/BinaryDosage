#' Function to read in SNPs from a binary dosage file
#'
#' Function to read in SNPs from a binary dosage file
#'
#' @param fileInfo
#' Information about the binary dosage file return from
#' GetBDoseInfo
#' @param SNPs
#' Vector of either SNP names or indices in binary dosage
#' file to extract values for
#' @param Subjects
#' Vector of either Subject IDs of indices in binary dosage
#' file to extract values for
#' @param geneProb
#' Include genetic probabilities - default true
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetSNPValues <- function(fileInfo, SNPs, Subjects, geneProb = TRUE) {
  if (missing(fileInfo))
    return (NULL)
  if (missing(SNPs)) {
    return (NULL)
  } else {
    if (is.vector(SNPs, "character")) {
      snpVec <- match(SNPs, fileInfo$SNPs$SNPID)
    } else if (is.vector(SNPs, "integer")) {
      snpVec <- SNPs
    } else {
      return (NULL)
    }
  }
  if (missing(Subjects)) {
    subVec <- 1L:nrow(fileInfo$Samples)
  } else {
    if (is.vector(Subjects, "character")) {
      subVec <- match(Subjects, fileInfo$Samples$SID)
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

  if (fileInfo$version == 1)
    geneProb <- FALSE
  if(geneProb == TRUE) {
    valueMatrix <- matrix(0, length(subVec), 4 * length(snpVec))
    gprob = 1
  } else {
    valueMatrix <- matrix(0, length(subVec), length(snpVec))
    gprob = 0
  }

  if (GetSNPValuesC(fileInfo$filename, fileInfo$filetype, gprob, subVec, snpVec,
                    fileInfo$Indices, valueMatrix, fileInfo$NumSamples, fileInfo$NumSNPs) == 1)
    return (NULL)

  if (geneProb == TRUE)
    colnames(valueMatrix) <- paste0(rep(fileInfo$SNPs$SNPID[snpVec], each = 4), rep(c("_dose", "_p0", "_p1", "_p2"), length(snpVec)))
  else
    colnames(valueMatrix) <- fileInfo$SNPs$SNPID[snpVec]
  rownames(valueMatrix) <- fileInfo$Samples$SID[subVec]

  return (valueMatrix)
}
