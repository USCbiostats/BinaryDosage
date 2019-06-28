AllocateBinaryDosageWriteMemory <- function(format, subformat, filename, funcData) {
  allocateFunc <- list(f1 <- c(AllocateBinaryDosageWriteMemory11, AllocateBinaryDosageWriteMemory12),
                       f2 <- c(AllocateBinaryDosageWriteMemory21, AllocateBinaryDosageWriteMemory22),
                       f3 <- c(AllocateBinaryDosageWriteMemory31, AllocateBinaryDosageWriteMemory32, AllocateBinaryDosageWriteMemory33, AllocateBinaryDosageWriteMemory34),
                       f4 <- c(AllocateBinaryDosageWriteMemory41, AllocateBinaryDosageWriteMemory42, AllocateBinaryDosageWriteMemory43, AllocateBinaryDosageWriteMemory44))
  return (allocateFunc[[format]][[subformat]](filename, funcData))
}

AllocateBinaryDosageWriteMemory11 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSubjects),
              usdosage = integer(funcData$numSubjects)))
}

AllocateBinaryDosageWriteMemory12 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSubjects),
              p0 = numeric(funcData$numSubjects),
              p1 = numeric(funcData$numSubjects),
              p2 = numeric(funcData$numSubjects),
              usdosage = integer(funcData$numSubjects),
              usp0 = integer(funcData$numSubjects),
              usp1 = integer(funcData$numSubjects),
              usp2 = integer(funcData$numSubjects)))
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
              dosage = numeric(funcData$numSubjects),
              usdosage = integer(funcData$numSubjects),
              altAlleleFreq = numeric(funcData$numSNPs)))
}

AllocateBinaryDosageWriteMemory42 <- function(filename, funcData) {
  return(list(filename = filename,
              dosage = numeric(funcData$numSubjects),
              p0 = numeric(funcData$numSubjects),
              p1 = numeric(funcData$numSubjects),
              p2 = numeric(funcData$numSubjects),
              usdosage = integer(funcData$numSubjects),
              usp0 = integer(funcData$numSubjects),
              usp1 = integer(funcData$numSubjects),
              usp2 = integer(funcData$numSubjects),
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
  return (WriteBinaryDosageP1Data(funcData$filename, funcData$dosage, funcData$p1, funcData$usdosage, funcData$usp1, 0))
}

WriteBinaryDosageData21 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 1))
}

WriteBinaryDosageData22 <- function(funcData) {
  return (WriteBinaryDosageP1Data(funcData$filename, funcData$dosage, funcData$p1, funcData$usdosage, funcData$usp1, 1))
}

WriteBinaryDosageData31 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 1))
}

WriteBinaryDosageData32 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData33 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 1))
}

WriteBinaryDosageData34 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData41 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 1))
}

WriteBinaryDosageData42 <- function(funcData) {
  return (0)
}

WriteBinaryDosageData43 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename, funcData$dosage, funcData$usdosage, 1))
}

WriteBinaryDosageData44 <- function(funcData) {
  return (0)
}

