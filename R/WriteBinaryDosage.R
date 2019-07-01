WriteFamilyAndMapFiles <- function(filename, samples, snps) {
  saveRDS(samples, filename[2])
  saveRDS(snps, filename[3])
  return (md5 <- c(digest(samples, "md5"), digest(snps, "md5")))
}
# Writes the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
WriteBinaryDosageHeader <- function(format, subformat, filename, funcData) {
  writeHeaderFunc <- list(f1 <- c(WriteBinaryDosageHeader11, WriteBinaryDosageHeader12),
                          f2 <- c(WriteBinaryDosageHeader21, WriteBinaryDosageHeader22),
                          f3 <- c(WriteBinaryDosageHeader31, WriteBinaryDosageHeader32, WriteBinaryDosageHeader33, WriteBinaryDosageHeader34),
                          f4 <- c(WriteBinaryDosageHeader41, WriteBinaryDosageHeader42, WriteBinaryDosageHeader43, WriteBinaryDosageHeader44))
  WriteBinaryDosageBaseHeader(filename[1], format - 1, subformat - 1)
  return (writeHeaderFunc[[format]][[subformat]](filename, funcData))
}

WriteBinaryDosageHeader11 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (md5)
#  return (WriteBinaryDosageHeader11C(filename))
}

WriteBinaryDosageHeader12 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (md5)
#  return (WriteBinaryDosageHeader12C(filename))
}

WriteBinaryDosageHeader21 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (md5)
#  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader22 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (md5)
#  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader31 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3A(filename[1], funcData$numSamples)
  return (md5)
  #  return (WriteBinaryDosageHeader31C(filename))
}

WriteBinaryDosageHeader32 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3A(filename[1], funcData$numSamples)
  return (md5)
  #  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader33 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2])
  return (md5)
  #  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader34 <- function(filename, funcData) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2])
  return (md5)
  #  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader41 <- function(filename, funcData) {
  WriteBinaryDosageHeader4A(filename[1], funcData$numSamples, funcData$numSNPs)
  WriteBDGroups(filename[1], funcData$numSamples)
  return (0)
  #  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader42 <- function(filename, funcData) {
  WriteBinaryDosageHeader4A(filename[1], funcData$numSamples, funcData$numSNPs)
  return (0)
  #  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader43 <- function(filename, funcData) {
  WriteBinaryDosageHeader4B(filename[1])
  return (0)
  #  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader44 <- function(filename, funcData) {
  WriteBinaryDosageHeader4B(filename[1])
  return (0)
  #  return (WriteBinaryDosageHeader22C(filename))
}
# Allocates memory needed to write binary dosage files
# Many of the routines do the same thing. Those that are copies
# just call the routine they are a copy of. This makes the code
# a little easier to read.
AllocateBinaryDosageWriteMemory <- function(format, subformat, filename, funcData) {
  allocateFunc <- list(f1 <- c(AllocateBinaryDosageWriteMemory11, AllocateBinaryDosageWriteMemory12),
                       f2 <- c(AllocateBinaryDosageWriteMemory21, AllocateBinaryDosageWriteMemory22),
                       f3 <- c(AllocateBinaryDosageWriteMemory31, AllocateBinaryDosageWriteMemory32, AllocateBinaryDosageWriteMemory33, AllocateBinaryDosageWriteMemory34),
                       f4 <- c(AllocateBinaryDosageWriteMemory41, AllocateBinaryDosageWriteMemory42, AllocateBinaryDosageWriteMemory43, AllocateBinaryDosageWriteMemory44))
  return (allocateFunc[[format]][[subformat]](filename, funcData))
}

AllocateBinaryDosageWriteMemory11 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSamples),
              usdosage = integer(funcData$numSamples)))
}

AllocateBinaryDosageWriteMemory12 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSamples),
              p0 = numeric(funcData$numSamples),
              p1 = numeric(funcData$numSamples),
              p2 = numeric(funcData$numSamples),
              usdosage = integer(funcData$numSamples),
              usp0 = integer(funcData$numSamples),
              usp1 = integer(funcData$numSamples),
              usp2 = integer(funcData$numSamples)))
}

AllocateBinaryDosageWriteMemory21 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory11(filename, funcData)
}

AllocateBinaryDosageWriteMemory22 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory12(filename, funcData)
}

AllocateBinaryDosageWriteMemory31 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory11(filename, funcData)
}

AllocateBinaryDosageWriteMemory32 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory12(filename, funcData)
}

AllocateBinaryDosageWriteMemory33 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory11(filename, funcData)
}

AllocateBinaryDosageWriteMemory34 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory12(filename, funcData)
}

AllocateBinaryDosageWriteMemory41 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSamples),
              usdosage = integer(funcData$numSamples),
              altAlleleFreq = numeric(funcData$numSNPs)))
}

AllocateBinaryDosageWriteMemory42 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSamples),
              p0 = numeric(funcData$numSamples),
              p1 = numeric(funcData$numSamples),
              p2 = numeric(funcData$numSamples),
              usdosage = integer(funcData$numSamples),
              usp0 = integer(funcData$numSamples),
              usp1 = integer(funcData$numSamples),
              usp2 = integer(funcData$numSamples),
              altAlleleFreq = numeric(funcData$numSNPs)))
}

AllocateBinaryDosageWriteMemory43 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory41(filename, funcData)
}

AllocateBinaryDosageWriteMemory44 <- function(filename, funcData) {
  AllocateBinaryDosageWriteMemory42(filename, funcData)
}

# Write binary dosage data at the end of the file
# Header has already been written
# funcData was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageFileData <- function(format, subformat, funcData) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageData11, WriteBinaryDosageData12),
                    f2 <- c(WriteBinaryDosageData21, WriteBinaryDosageData22),
                    f3 <- c(WriteBinaryDosageData31, WriteBinaryDosageData32, WriteBinaryDosageData33, WriteBinaryDosageData34),
                    f4 <- c(WriteBinaryDosageData41, WriteBinaryDosageData42, WriteBinaryDosageData43, WriteBinaryDosageData44))
  return (writeFunc[[format]][[subformat]](funcData))
}

WriteBinaryDosageData11 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 0))
}

WriteBinaryDosageData12 <- function(funcData) {
  return (WriteBinaryDosageP1Data(funcData$filename, funcData$dosage, funcData$p1, funcData$usdosage, funcData$usp1, 1))
}

WriteBinaryDosageData21 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData22 <- function(funcData) {
  return (WriteBinaryDosageP1Data(funcData$filename, funcData$dosage, funcData$p1, funcData$usdosage, funcData$usp1, 2))
}

WriteBinaryDosageData31 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData32 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData33 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData34 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData41 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData42 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData43 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData44 <- function(funcData) {
  return (0)
}

