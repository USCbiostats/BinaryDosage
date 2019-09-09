bdapply <- function(bdinfo, func, ...) {
  retval <- vector("list", nrow(bdinfo$snps))
  dosage <- numeric(nrow(bdinfo$samples))
  p0 <- numeric(nrow(bdinfo$samples))
  p1 <- numeric(nrow(bdinfo$samples))
  p2 <- numeric(nrow(bdinfo$samples))
  us <- integer(2 * nrow(bdinfo$samples))
  for (i in 1:nrow(bdinfo$snps)) {
    dosage[1:nrow(bdinfo$samples)] <- NA
    p0[1:nrow(bdinfo$samples)] <- NA
    p1[1:nrow(bdinfo$samples)] <- NA
    p2[1:nrow(bdinfo$samples)] <- NA
    ReadBinaryDosageData(bdinfo, i, dosage, p0, p1, p2, us)
    retval[[i]] <- func(dosage, p0, p1, p2, ...)
  }
  return (retval)
}
