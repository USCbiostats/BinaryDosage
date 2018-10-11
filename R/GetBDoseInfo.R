#' Function to read in information about binary dosage file
#'
#' Function to read in information about binary dosage file.
#' This information is used when read the dosage and genetic
#' probabilities from the file
#'
#' @param bdFile
#' Name of the binary dosage file
#' @param famFile
#' Name of file with subject data. This file is in the format
#' used by plink .fam files. Not needed for format 4 or greater
#' binary dosage files.
#' @param mapFile
#' Name of file with SNP information. This file is in the format
#' used by plink for extended map file files. Not needed for
#' format 4 or greater binary dosage files.
#' @param index
#' Index the SNPs for faster reading.
#' @param sep
#' Separator in fam and map files
#' @return
#' List with information about the file including subject and
#' SNP information
#' @export
GetBDoseInfo <- function(bdFile, famFile = "", mapFile = "", index = FALSE, sep = '\t') {
  if (missing(bdFile))
    return (NULL)
  if (index == TRUE)
    x = 1L
  else
    x = 0L

  bdf <- GetBDoseFormatC(bdFile)
  if (bdf$Format == 4)
    return (GetBinaryDosage4Info(bdFile, x))

  subjects <- read.table(famFile, sep=sep, stringsAsFactors = FALSE, nrows = 1)
  if (ncol(subjects) < 5 || ncol(subjects) > 6) {
    print("Fam file does not have correct number of columns")
    return (list())
  }
  if (ncol(subjects) == 5) {
    subjects <- as.data.frame(read.table(famFile, sep = sep, stringsAsFactors = FALSE,
                                         colClasses = c("character","NULL", "NULL","NULL", "NULL")))
    subjects$FID = ""
    subjects <- subjects[,c(2,1)]
  } else {
    subjects <- read.table(famFile, sep = sep, stringsAsFactors = FALSE,
                           colClasses = c("character", "character", "NULL", "NULL","NULL", "NULL"))
  }
  colnames(subjects) <- c("FID", "SID")
  snps <- read.table(mapFile, sep=sep, stringsAsFactors = FALSE, nrows = 1)
  if (ncol(snps) != 6) {
    print("Map file dose not have 6 columns")
    return (list())
  }
  snps <- read.table(mapFile, sep=sep, stringsAsFactors = FALSE,
                     colClasses = c("character", "character", "NULL", "integer","character", "character"))
  colnames(snps) <- c("CHR", "SNPID", "LOC", "REF", "ALT")
  return (GetBinaryDosage1Info(bdFile, subjects, snps, x))
}
