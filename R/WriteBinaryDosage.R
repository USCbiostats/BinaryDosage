WriteBDTest <- function(format, subformat, funcData) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageData11, WriteBinaryDosageData12),
                    f2 <- c(WriteBinaryDosageData21, WriteBinaryDosageData22),
                    f3 <- c(WriteBinaryDosageData31, WriteBinaryDosageData33),
                    f4 <- c(WriteBinaryDosageData41, WriteBinaryDosageData43))
  return (writeFunc[[format]][[subformat]](funcData))
}

WriteBinaryDosageData11 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                0))
}

WriteBinaryDosageData12 <- function(funcData) {
  return (WriteBinaryDosageP1Data(funcData$filename,
                                  funcData$dosage,
                                  funcData$p1,
                                  funcData$usdosage,
                                  funcData$usp1,
                                  0))
}

WriteBinaryDosageData21 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                1))
}

WriteBinaryDosageData22 <- function(funcData) {
  return (WriteBinaryDosageP1Data(funcData$filename,
                                  funcData$dosage,
                                  funcData$p1,
                                  funcData$usdosage,
                                  funcData$usp1,
                                  1))
}

WriteBinaryDosageData31 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                1))
}

WriteBinaryDosageData33 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                1))
}

WriteBinaryDosageData41 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                1))
}

WriteBinaryDosageData43 <- function(funcData) {
  return (WriteBinaryDosageData(funcData$filename,
                                funcData$dosage,
                                funcData$usdosage,
                                1))
}
