bdapply <- function(bdinfo, func, funcdata) {
  dosage <- numeric(bdinfo$numSamples)
  p0 <- numeric(bdinfo$numSamples)
  p1 <- numeric(bdinfo$numSamples)
  p2 <- numeric(bdinfo$numSamples)
  us <- integer(2*bdinfo$numSamples)
  for (i in 1:bdinfo$numSNPs) {
    ReadBinaryDosageData(bdinfo, i, dosage, p0, p1, p2, us)
    func(funcdata, dosage, p0, p1, p2)
  }
}
