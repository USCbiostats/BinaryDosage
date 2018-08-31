#' Function to merge binary dosage files in format 4.2
#'
#' Function to merge binary dosage files in format 4.2
#'
#' @param mergeFile
#' Name of file that will contain the merged data.
#' @param filesToMerge
#' A vector of strings. Names of the files to merge.
#' @return
#' 0 - Successfully merged
#' 1 - Merge failed
#' @export
MergeBD42 <- function(mergedFile, filesToMerge) {
  if (length(filesToMerge) < 2) {
    print("At least two files must be specified to merge")
    return (1)
  }
  fid <- character();
  iid <- character();
  matchList <- vector("list", length = length(filesToMerge))
  for (i in c(1:length(filesToMerge))) {
    bdInfo <- GetBinaryDosageInfo(filesToMerge[i])
    if (bdInfo$format != 4 || bdInfo$version != 2) {
      print ("Not all files in format 4.2")
      return (1)
    }
    fid <- c(fid, bdInfo$Samples$FID)
    iid <- c(iid, bdInfo$Samples$SID)
    if (i == 1)
      snpList1 <- bdInfo$SNPs
    snpList2 <- bdInfo$SNPs
    matchList[[i]] = row.match(snpList1, snpList2)
  }
  subjects <- data.frame(fid, iid, stringsAsFactors = FALSE);
  if (length(unique(subjects)) != length(subjects)) {
    print("Sample IDs are not unique")
    return (list())
  }
  snpLocations <- data.frame(matchList)
  snpLocations <- as.matrix(snpLocations[complete.cases(snpLocations),])
  snpsToMerge <- snpList1[snpLocations[,1],]
  if (nrow(snpsToMerge) == 0) {
    print("Intersection of SNPs is empty")
    return (list())
  }

  mergeInfo <- list(subjects = subjects, locations = snpLocations, snpsToMerge = snpsToMerge)

  return (Merge42C(mergedFile, filesToMerge, mergeInfo))
}
