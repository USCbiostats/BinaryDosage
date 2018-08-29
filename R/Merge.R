GetSubjects <- function(bdInfoList) {
  fid <- character();
  iid <- character();
  for (i in c(1:(length(bdInfoList)))) {
    fid <- c(fid, bdInfoList[[i]]$Samples$FID)
    iid <- c(iid, bdInfoList[[i]]$Samples$SID)
  }
  df <- data.frame(fid, iid, stringsAsFactors = FALSE);
  if (length(unique(df)) != length(df))
    return (data.frame())
  return (df)
}

CombineSubjects <- function(bdInfoList) {
  subID <- numeric()
  for (i in c(1:(length(bdInfoList))))
    subID <- c(subID, bdInfoList[[i]]$subjects$IID)
  return (subID)
}

SNPIntersection <- function(bdInfoList, snpsToUse) {
  df <- bdInfoList[[1]]$SNPs[snpsToUse[,1],]
  snpID <- paste(df$CHR, df$BP, sep = ":")
  noSNPID = all(df$SNP == snpID)
  singleChromosome = (length(unique(df$CHR)) == 1)
  return (list(SNPs = df, noSNPID = noSNPID, singleChromosome = singleChromosome))
}

GetSNPList <- function(bdInfoList) {
  matchList <- vector("list", length = length(bdInfoList))
  snpList1 <- paste(bdInfoList[[1]]$SNPs$SNP, bdInfoList[[1]]$SNPs$A1, bdInfoList[[1]]$SNPs$A2)
  for (i in c(1:length(bdInfoList))) {
    snpList2 <- paste(bdInfoList[[i]]$SNPs$SNP, bdInfoList[[i]]$SNPs$A1, bdInfoList[[i]]$SNPs$A2)
    matchList[[i]] <- match(snpList1, snpList2)
  }
  cn <- paste("set", 1:length(bdInfoList), sep = '_')
  df <- data.frame(matchList)
  colnames(df) <- cn
  df <- df[complete.cases(df),]
  return (df)
}

MergeBD <- function(filenames) {
  bdInfoList <- vector("list", length = length(filenames))
  for (i in c(1:length(filenames)))
    bdInfoList[[i]] <- GetBinaryDosageInfo(filenames[i])
  return (bdInfoList)
}
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
  bdInfoList <- vector("list", length = length(filesToMerge))
  for (i in c(1:length(filesToMerge))) {
    bdInfoList[[i]] <- GetBinaryDosageInfo(filesToMerge[i])
    if (bdInfoList[[i]]$format != 4 || bdInfoList[[i]]$version != 2) {
      print ("Not all files in format 4.2")
      return (1)
    }
  }
  subjects <- GetSubjects(bdInfoList = bdInfoList)
  if (nrow(subjects) == 0) {
    print ("Duplicate subject IDs exist")
    return (2)
  }
  snpLocations <- GetSNPList(bdInfoList = bdInfoList)
  snpsToMerge <- SNPIntersection(bdInfoList = bdInfoList, snpsToUse = snpLocations)
  if (nrow(snpsToMerge$SNPs) == 0) {
    print ("Intersection of SNPs has size 0")
    return (3)
  }

  snpName <- paste(snpsToMerge$SNPs$Chromosome, snpsToMerge$SNPs$Location, sep = ':')
  if (all(snpName == snpsToMerge$SNPs$SNPID))
    snpsToMerge$SNPs$SNPID <- ""

  return (list(bdInfo = bdInfoList, subjects = subjects, locations = as.matrix(snpLocations), snpsToMerge = snpsToMerge))
}
