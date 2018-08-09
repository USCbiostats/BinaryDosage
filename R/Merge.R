RepeatSubjects <- function(bdInfoList) {
  for (i in c(1:(length(bdInfoList)- 1))) {
    fs1 <- paste(bdInfoList[[i]]$subjects$FID, bdInfoList[[i]]$subjects$IID, sep = '_')
    for (j in c((i + 1):length(bdInfoList))) {
      fs2 <- paste(bdInfoList[[j]]$subjects$FID, bdInfoList[[j]]$subjects$IID, sep = '_')
      if(any(match(fs1, fs2), na.rm = TRUE) == TRUE)
        return (TRUE)
    }
  }
  return(FALSE)
}

CombineSubjects <- function(bdInfoList) {
  subID <- numeric()
  for (i in c(1:(length(bdInfoList))))
    subID <- c(subID, bdInfoList[[i]]$subjects$IID)
  return (subID)
}

SNPIntersection <- function(bdInfoList, snpsToUse) {
  df <- bdInfoList[[1]]$snps[snpsToUse[,1],]
  snpID <- paste(df$CHR, df$BP, sep = ":")
  noSNPID = all(df$SNP == snpID)
  singleChromosme = (length(unique(df$CHR)) == 1)
  return (list(df, noSNPID, singleChromosme))
}

GetSNPList <- function(bdInfoList) {
  matchList <- vector("list", length = length(bdInfoList))
  snpList1 <- paste(bdInfoList[[1]]$snps$SNP, bdInfoList[[1]]$snps$A1, bdInfoList[[1]]$snps$A2)
  for (i in c(1:length(bdInfoList))) {
    snpList2 <- paste(bdInfoList[[i]]$snps$SNP, bdInfoList[[i]]$snps$A1, bdInfoList[[i]]$snps$A2)
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
