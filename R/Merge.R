#' Function to merge binary dosage files
#'
#' Function to merge binary dosage files
#'
#' @param mergeFile
#' Name of file that will contain the merged data.
#' @param filesToMerge
#' A vector of strings. Names of the files to merge.
#' @param famFilesToMerge
#' A vector of strings. These are the names of the fam files
#' associated with files to merge. If the files to merge are
#' in format 4 or later, this is not needed.
#' @param mapFilesToMerge
#' A vector of strings. These are the names of the map files
#' associated with files to merge. If the files to merge are
#' in format 4 or later, this is not needed.
#' @param mergedFamFile
#' Name of fam file to be associated with the binary dosage
#' file. This is only needed if format is 3 or less.
#' @param mergedMapFile
#' Name of map file to be associated with the binary dosage
#' file. This is only needed if format is 3 or less.
#' @param format
#' Format of the merged binary dosage file. Default value is 4.
#' It is recommended to use the default value.
#' @param version
#' Version of the format to use. Default value is 2. If the files
#' being merged do not contain gene probabilities, change this
#' value to 1.
#' @param batchSize
#' Number of SNPs to process at a time
#' @return
#' 0 - Successfully merged
#' 1 - Merge failed
#' @export
MergeBD <- function(mergedFile, filesToMerge, famFilesToMerge, mapFilesToMerge,
                      mergedFamFile = "", mergeMapFile = "", format = 4, version = 2, batchSize = 1000) {
  print("Entering")
  if (length(filesToMerge) < 2) {
    print("At least two files must be specified to merge")
    return (1)
  }
  if (missing(famFilesToMerge) == FALSE) {
    if (missing(mapFilesToMerge) == FALSE) {
      if (length(famFilesToMerge) != length(filesToMerge) | length(mapFilesToMerge) != length(filesToMerge)) {
        print("Lengths of filesToMerge, famFilesToMerge, and mapFilesToMerge are not the same")
        return (1)
      }
    } else {
      print("If either famFilesToMerge or mapFilesToMerge is provided, the other must be provided")
      print(1)
    }
  } else if (missing(mapFilesToMerge) == FALSE) {
    print("If either famFilesToMerge or mapFilesToMerge is provided, the other must be provided")
    print(1)
  }
  if (format < 1 | format > 4 | version < 1 | version > 2) {
    print("Format and version are not valid")
    return (1)
  }

  print("Before fid")
  fid <- character();
  iid <- character();
  matchList <- vector("list", length = length(filesToMerge))
  bdInfoList <- vector("list", length = length(filesToMerge))

  for (i in c(1:length(filesToMerge))) {
    print("Before bdInfo")
    if (missing(famFilesToMerge) == TRUE) {
      bdInfo <- GetBDoseInfo(filesToMerge[i], index = 1)
    } else {
      bdInfo <- GetBDoseInfo(filesToMerge[i], famFilesToMerge[i], mapFilesToMerge[i])
    }
    print("After bdInfo")
    if (bdInfo$version == 1 & version == 2) {
      print ("Cannot merge files in version 1 to a file with version 2")
      return (1)
    }
    fid <- c(fid, bdInfo$Samples$FID)
    iid <- c(iid, bdInfo$Samples$SID)
    if (i == 1)
      snpList1 <- bdInfo$SNPs
    snpList2 <- bdInfo$SNPs
    matchList[[i]] = row.match(snpList1, snpList2)
    bdInfoList[[i]] = bdInfo
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

  mergeInfo <- list(subjects = subjects, locations = snpLocations, snpsToMerge = snpsToMerge, bdInfoList = bdInfoList)
  return (MergeBDC(mergedFile, filesToMerge, mergeInfo, bdInfoList, mergedFamFile, mergeMapFile, format, version, batchSize))
#  return (mergeInfo)
}
