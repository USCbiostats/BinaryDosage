ReadVCFSection <- function(vcfInfo, firstSNP, snpValues) {
  if (vcfInfo$gzipped == FALSE) {
    con <- file(vcfInfo$filename, "r")
    seek(con, sum(vcfInfo$Indices[1:firstSNP]))
  } else {
    con <- gzfile(vcfInfo$filename, "r")
    line <- readLines(con, n = vcfInfo$headersize)
    if (firstSNP > 1)
      line <- readLines(con, n = firstSNP - 1)
  }

  line <- readLines(con, n = 1)
  x <- unlist(strsplit(line, "\t"))
  y <- unlist(strsplit(x[10:length(x)], ":"))
  if (length(vcfInfo$additionalInfo$dataColumns$dosage) == 1) {
    dosageCol <- vcfInfo$additionalInfo$dataColumns$dosage
    gpCol <- vcfInfo$additionalInfo$dataColumns$genotypeProb
    numValues <- vcfInfo$additionalInfo$dataColumns$numValues
  } else {
    dosageCol <- vcfInfo$additionalInfo$dataColumns$dosage[firstSNP]
    gpCol <- vcfInfo$additionalInfo$dataColumns$genotypeProb[firstSNP]
    numValues <- vcfInfo$additionalInfo$dataColumns$numValues[firstSNP]
  }
  if(is.na(dosageCol) == FALSE) {
    dosage <- as.numeric(y[seq(dosageCol, length(y) - numValues + dosageCol, numValues)])
  }
  if(is.na(gpCol) == FALSE) {
    gpString <- y[seq(gpCol, length(y) - numValues + gpCol, numValues)]
    z <- unlist(strsplit(gpString, ","))
    p0 <- as.numeric(z[seq(1, length(z) - 2, 3)])
    p1 <- as.numeric(z[seq(2, length(z) - 1, 3)])
    p2 <- as.numeric(z[seq(3, length(z), 3)])
  }
  close(con)

  return (data.frame(dosage, p0, p1, p2, stringsAsFactors = FALSE))
}
