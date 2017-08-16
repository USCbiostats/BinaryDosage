#' @importFrom Rcpp evalCpp
#' @importFrom data.table fread setkey data.table
#' @useDynLib GxEScanR
NULL

#' Function to read in information about binary dosage file
#'
#' Function to read in information about binary dosage file.
#' This information is used when read the dosage and genetic
#' probabilities from the file
#'
#' @param filename
#' Name of the binary dosage file
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetBinaryFileInfo <- function(filename) {
  return (ReadBDInfo(filename))
}

#' Function to create data table to contain dosage values
#'
#' Function to create data table to contain dosage values
#' This data table is passed to the dosage extraction routines.
#' The dosage values are then added to this data table
#'
#' @param binaryFileInfo
#' A list returned from GetBinaryFileInfo
#' @return
#' A data.table with the subject the information
#' @export
CreateSubjectDataTable <- function(binaryFileInfo) {
  return (data.table(binaryFileInfo$Subjects))
}

#' Function to read SNPs from a binary dosage file
#'
#' Function to read SNPs from a binary dosage file
#' GetBinaryFileInfo and CreateSubjectDataTable need
#' to be run first
#'
#' @param binaryFileInfo
#' A list returned from GetBinaryFileInfo
#' @param subjectData
#' Data table with subject data - returned from CreateSubjectDataTable
#' @param SNPNumbers
#' Integer vector on SNPs to retrieve
#' @param maxSNPs
#' Maximum number of SNPs to read in
#' @return
#' List with information about what SNPs were read in
#' This list may be pass to ReadMoreSNPs to retrieve additional
#' SNPs if the number of SNPs specified in SNPNumbers is larger
#' that maxSNPs
#' @export
ReadBinarySNPData <- function(binaryFileInfo, subjectData, SNPNumbers, maxSNPs) {
  if (length(SNPNumbers) < maxSNPs)
    maxSNPs <- length(SNPNumbers)
  firstSNP <- 1
  sc <- firstSNP:maxSNPs
  d <- rep("Dose", maxSNPs)
  d <- paste(d, as.character(sc), sep = "_")
  p0 <- rep("p0", maxSNPs)
  p0 <- paste(p0, as.character(sc), sep = "_")
  p1 <- rep("p1", maxSNPs)
  p1 <- paste(p1, as.character(sc), sep = "_")
  p2 <- rep("p2", maxSNPs)
  p2 <- paste(p2, as.character(sc), sep = "_")
  subjectData[, c(d) := NA_real_]
  subjectData[, c(p0) := NA_real_]
  subjectData[, c(p1) := NA_real_]
  subjectData[, c(p2) := NA_real_]
  ptrd <- rep(as.numeric(0), maxSNPs)
  ptrp0 <- rep(as.numeric(0), maxSNPs)
  ptrp1 <- rep(as.numeric(0), maxSNPs)
  ptrp2 <- rep(as.numeric(0), maxSNPs)
  for (i in 1:maxSNPs) {
    eval(parse(text = paste("FindPointer(subjectData$",d[i], ", ptrd,", as.character(i), ")", sep = "")))
    eval(parse(text = paste("FindPointer(subjectData$",p0[i], ", ptrp0,", as.character(i), ")", sep = "")))
    eval(parse(text = paste("FindPointer(subjectData$",p1[i], ", ptrp1,", as.character(i), ")", sep = "")))
    eval(parse(text = paste("FindPointer(subjectData$",p2[i], ", ptrp2,", as.character(i), ")", sep = "")))
  }
#  PrintPointer(ptrd)
#  PrintPointer(ptrp0)
#  PrintPointer(ptrp1)
#  PrintPointer(ptrp2)
  ReadBDSNPs(binaryFileInfo$FileData$FileName, binaryFileInfo$FileData$NumSubject, binaryFileInfo$FileData$NumSNPs,
             binaryFileInfo$FileData$DosageStart, binaryFileInfo$FileData$DataSize, SNPNumbers[sc], ptrd, ptrp0, ptrp1, ptrp2)
  return (list(firstSNP = firstSNP, lastSNP = maxSNPs,
               ptrd = ptrd, ptrp0 = ptrp0, ptrp1 = ptrp1, ptrp2 = ptrp2,
               SNPNumbers = SNPNumbers, maxSNPs = maxSNPs, moreSNPs = (maxSNPs != length(SNPNumbers))))
}

#' Function to read SNPs from a binary dosage file
#'
#' Function to read SNPs from a binary dosage file
#' GetBinaryFileInfo and CreateSubjectDataTable need
#' to be run first
#'
#' @param binaryFileInfo
#' Information about binary file returned from GetBinaryFileInfo
#' @param readInfo
#' A list returned from ReadBinarySNPData
#' @return
#' List with information about what SNPs were read in
#' This list may be pass to ReadMoreBinarySNPData to retrieve additional
#' SNPs if the number of SNPs specified in SNPNumbers is larger
#' that maxSNPs
#' @export
ReadMoreBinarySNPData <- function(binaryFileInfo, readInfo) {
  if (readInfo$moreSNPs == FALSE)
    return (readInfo)

  readInfo$firstSNP <- readInfo$lastSNP + 1
  if (length(readInfo$SNPNumbers) < readInfo$firstSNP + readInfo$maxSNPs) {
    readInfo$lastSNP <- length(readInfo$SNPNumbers)
    readInfo$moreSNPs <- FALSE
  } else {
    readInfo$lastSNP <- readInfo$firstSNP + readInfo$maxSNPs - 1
  }
#  print(length(readInfo$SNPNumbers))
#  print(firstSNP)
#  print(readInfo$lastSNP)
  sc <- readInfo$firstSNP:readInfo$lastSNP
  ReadBDSNPs(binaryFileInfo$FileData$FileName, binaryFileInfo$FileData$NumSubject, binaryFileInfo$FileData$NumSNPs,
             binaryFileInfo$FileData$DosageStart, binaryFileInfo$FileData$DataSize, readInfo$SNPNumbers[sc],
             readInfo$ptrd, readInfo$ptrp0, readInfo$ptrp1, readInfo$ptrp2)
  return (readInfo)
}

#' Function to loop through the SNPs in a binary dosage file
#'
#' Function to loop through the SNPs in a binary dosage file
#' and apply a function to each each SNP
#'
#' @param binaryDosageInfo
#' Information returned from GetBinaryFileInfo
#' @param function2apply
#' Function to be preformed on each SNP
#' @param outputValues
#' Place to store the results of the function
#' This is usually a data.table because it can be modified by
#' function
#' @param extraData
#' Addition data to be passed to the function to apply
#' @param snps2use
#' List of indices of the SNPs in the binary dosage file to use
#' @param maxSNPs
#' Maximum number of SNPs to read in at one time - current limit is 200
#' This needs work. It is an issue with data.tables
#' @return
#' 0 - success - currently always succeeds
#' 1 - failure
#' The results of the function to apply are stored in outputValues
#' @export
SNPApply <- function(binaryDosageInfo, function2apply, outputValues, extraData = NULL, snps2use, maxSNPs) {
  fileSubjects <- CreateSubjectDataTable(binaryDosageInfo)
  SNPStatus <- ReadBinarySNPData(binaryDosageInfo, fileSubjects, snps2use, maxSNPs = maxSNPs)
  firstCol <- which(names(fileSubjects) %in% "Dose_1")
  lastCol <- which(names(fileSubjects) %in% paste("Dose", as.character(SNPStatus$lastSNP - SNPStatus$firstSNP + 1), sep = "_"))
  function2apply(fileSubjects, outputValues, SNPStatus$firstSNP:SNPStatus$lastSNP, firstCol, lastCol, extraData)
  while (SNPStatus$moreSNPs) {
    SNPStatus <- ReadMoreBinarySNPData(binaryDosageInfo, SNPStatus)
    if (SNPStatus$lastSNP - SNPStatus$firstSNP + 1 != maxSNPs)
      lastCol <- which(names(fileSubjects) %in% paste("Dose", as.character(SNPStatus$lastSNP - SNPStatus$firstSNP + 1), sep = "_"))
    function2apply(fileSubjects, snpFreq, SNPStatus$firstSNP:SNPStatus$lastSNP, firstCol, lastCol, extraData)
  }
  return (0)
}
