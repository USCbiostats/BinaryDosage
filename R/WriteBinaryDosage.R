# Formats 1, 2, and 3 all have a separate family and map file
# This routine saves the data frames in RDS format
WriteFamilyAndMapFiles <- function(filename, samples, snps) {
  saveRDS(samples, filename[2])
  saveRDS(snps, filename[3])
  return (md5 <- c(digest(samples, "md5"), digest(snps, "md5")))
}

# Format 4 has the family data in the binary dosage file
WriteBDFamilyInfo <- function(filename, funcData, numsuboffset, suboffset, snpoffset) {
  sid <- paste(funcData$samples$SID, collapse = '\t')
  if (funcData$usesFID == TRUE)
    fid <- paste(funcData$samples$FID, collapse = '\t')
  else
    fid <- ""
  return (WriteBDFamilyInfoC(filename[1],
                             funcData$numSamples,
                             sid[1],
                             fid[1],
                             numsuboffset,
                             suboffset,
                             snpoffset))
}

FindSNPInfoString <- function(toFind, snpInfo, numSNPs) {
  infocol <- match(toFind, colnames(snpInfo))
  if (is.na(infocol) == FALSE)
    return (snpInfo[,infocol])
  return ("")
}

FindSNPInfoNumeric <- function(toFind, snpInfo, numSNPs, bdoptions) {
  if (is.na(match(toFind, bdoptions)) == FALSE)
    return (numeric(numSNPs))
  infocol <- match(toFind, colnames(snpInfo))
  if (is.na(infocol) == FALSE)
    return (snpInfo[,infocol])
  return (numeric(0))
}

# Format 4 has the SNP data in the binary dosage file
WriteBDSNPInfo <- function(filename, funcData, numsnpoffset, snpOptionsoffset, snpoffset, nextoffset, bdoptions) {
  snpcolnames <- colnames(funcData$SNPs)
  if (funcData$snpidformat == 0)
    snpid <- FindSNPInfoString("SNPID", funcData$SNPs, funcData$numSNPs)
  else
    snpid <- ""
  chr <- FindSNPInfoString("Chromosome", funcData$SNPs, funcData$numSNPs)
  if (funcData$onechr == TRUE)
    chr <- chr[1]
  if (is.na(match("Location", colnames(funcData$SNPs))) == TRUE)
    loc = rep(0L, funcData$numSNPs)
  else
    loc <- funcData$SNPs$Location
  ref <- FindSNPInfoString("Reference", funcData$SNPs, funcData$numSNPs)
  alt <- FindSNPInfoString("Alternate", funcData$SNPs, funcData$numSNPs)

  aaf <- FindSNPInfoNumeric("aaf", funcData$snpInfo, funcData$numSNPs, bdoptions)
  maf <- FindSNPInfoNumeric("maf", funcData$snpInfo, funcData$numSNPs, bdoptions)
  # Average call cannot be calculated, therefore bdoptions are meaningless
  avgCall <- FindSNPInfoNumeric("avgCall", funcData$snpInfo, funcData$numSNPs, "")
  rsq <- FindSNPInfoNumeric("rsq", funcData$snpInfo, funcData$numSNPs, bdoptions)

  snpid <- paste0(snpid, collapse = '\t')
  chr <- paste0(chr, collapse = '\t')
  ref <- paste0(ref, collapse = '\t')
  alt <- paste0(alt, collapse = '\t')
  WriteBDSNPInfoC(filename[1], funcData$numSNPs, snpid, chr, loc, ref, alt,
                  aaf, maf, avgCall, rsq, numsnpoffset, snpOptionsoffset, snpoffset, nextoffset)
  return (list(snpid = snpid,
               chr = chr,
               loc = loc,
               ref = ref,
               alt = alt,
               aaf = aaf,
               maf = maf,
               avgCall = avgCall,
               rsq = rsq))
}
# Writes the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
WriteBinaryDosageHeader <- function(format, subformat, filename, funcData, bdoptions) {
  writeHeaderFunc <- list(f1 <- c(WriteBinaryDosageHeader11, WriteBinaryDosageHeader12),
                          f2 <- c(WriteBinaryDosageHeader21, WriteBinaryDosageHeader22),
                          f3 <- c(WriteBinaryDosageHeader31, WriteBinaryDosageHeader32, WriteBinaryDosageHeader33, WriteBinaryDosageHeader34),
                          f4 <- c(WriteBinaryDosageHeader41, WriteBinaryDosageHeader42, WriteBinaryDosageHeader43, WriteBinaryDosageHeader44))
  WriteBinaryDosageBaseHeader(filename[1], format - 1, subformat - 1)
  return (writeHeaderFunc[[format]][[subformat]](filename, funcData, bdoptions))
}

WriteBinaryDosageHeader11 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (0)
#  return (WriteBinaryDosageHeader11C(filename))
}

WriteBinaryDosageHeader12 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (0)
#  return (WriteBinaryDosageHeader12C(filename))
}

WriteBinaryDosageHeader21 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (0)
#  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader22 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (0)
#  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader31 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3A(filename[1], funcData$numSamples)
  return (0)
  #  return (WriteBinaryDosageHeader31C(filename))
}

WriteBinaryDosageHeader32 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3A(filename[1], funcData$numSamples)
  return (0)
  #  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader33 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2])
  return (0)
  #  return (WriteBinaryDosageHeader21C(filename))
}

WriteBinaryDosageHeader34 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2])
  WriteBDIndexArray3_4(filename[1], funcData$numSNPs)
  return (0)
  #  return (WriteBinaryDosageHeader22C(filename))
}

WriteBinaryDosageHeader41 <- function(filename, funcData, bdoptions) {
  WriteBinaryDosageHeader4A(filename[1], funcData$numSamples, funcData$numSNPs)
  WriteBDGroups(filename[1], funcData$numSamples)
  WriteBDFamilyInfo(filename[1], funcData, 8, 28, 32)
  WriteBDSNPInfo(filename[1], funcData, 12, 24, 32, 36, bdoptions)
  return (0)
}

WriteBinaryDosageHeader42 <- function(filename, funcData, bdoptions) {
  return (WriteBinaryDosageHeader41(filename, funcData, bdoptions))
}

WriteBinaryDosageHeader4_3_4 <- function(filename, funcData, bdoptions) {
  WriteBinaryDosageHeader4B(filename[1], funcData$numSamples, funcData$numSNPs)
  WriteBDGroups2(filename[1], funcData$numSamples)
  WriteBDFamilyInfo(filename[1], funcData, -1, 8, 12)
  WriteBDSNPInfo(filename[1], funcData, -1, -1, 12, 16, bdoptions)
  return (0)
}

WriteBinaryDosageHeader43 <- function(filename, funcData, bdoptions) {
  WriteBinaryDosageHeader4_3_4(filename, funcData, bdoptions)
  return (0)
}

WriteBinaryDosageHeader44 <- function(filename, funcData, bdoptions) {
  WriteBinaryDosageHeader4_3_4(filename, funcData, bdoptions)
  WriteBDIndexArray4_4(filename[1], funcData$numSNPs, 16, 20)
  return (0)
}

# Allocates memory needed to write binary dosage files
# Many of the routines do the same thing. Those that are copies
# just call the routine they are a copy of. This makes the code
# a little easier to read.
AllocateBinaryDosageWriteMemory <- function(funcData) {
  allocateFunc <- list(f1 <- c(AllocateBinaryDosageWriteMemory11, AllocateBinaryDosageWriteMemory12),
                       f2 <- c(AllocateBinaryDosageWriteMemory21, AllocateBinaryDosageWriteMemory22),
                       f3 <- c(AllocateBinaryDosageWriteMemory31, AllocateBinaryDosageWriteMemory32, AllocateBinaryDosageWriteMemory33, AllocateBinaryDosageWriteMemory34),
                       f4 <- c(AllocateBinaryDosageWriteMemory41, AllocateBinaryDosageWriteMemory42, AllocateBinaryDosageWriteMemory43, AllocateBinaryDosageWriteMemory44))
  return (allocateFunc[[funcData$format]][[funcData$subformat]](funcData))
}

AllocateBinaryDosageWriteMemory11 <- function(funcData) {
  return(list(filename = funcData$filename,
              format = funcData$format,
              subformat = funcData$subformat,
              usdosage = integer(funcData$numSamples)))
}

AllocateBinaryDosageWriteMemory12 <- function(funcData) {
  return(list(filename = funcData$filename,
              format = funcData$format,
              subformat = funcData$subformat,
              usp1 = integer(funcData$numSamples),
              usp2 = integer(funcData$numSamples)))
}

AllocateBinaryDosageWriteMemory21 <- function(funcData) {
  AllocateBinaryDosageWriteMemory11(funcData)
}

AllocateBinaryDosageWriteMemory22 <- function(funcData) {
  AllocateBinaryDosageWriteMemory12(funcData)
}

AllocateBinaryDosageWriteMemory31 <- function(funcData) {
  AllocateBinaryDosageWriteMemory11(funcData)
}

AllocateBinaryDosageWriteMemory32 <- function(funcData) {
  AllocateBinaryDosageWriteMemory12(funcData)
}

AllocateBinaryDosageWriteMemory33 <- function(funcData) {
  AllocateBinaryDosageWriteMemory11(funcData)
}

AllocateBinaryDosageWriteMemory34 <- function(funcData) {
  AllocateBinaryDosageWriteMemory12(funcData)
}

AllocateBinaryDosageWriteMemory41 <- function(funcData) {
  return(list(filename = funcData$filename,
              format = funcData$format,
              subformat = funcData$subformat,
              usdosage = integer(funcData$numSamples),
              altAlleleFreq = numeric(funcData$numSNPs)))
}

AllocateBinaryDosageWriteMemory42 <- function(funcData) {
  return(list(filename = funcData$filename,
              format = funcData$format,
              subformat = funcData$subformat,
              usdosage = integer(funcData$numSamples),
              usp0 = integer(funcData$numSamples),
              usp1 = integer(funcData$numSamples),
              usp2 = integer(funcData$numSamples),
              altAlleleFreq = numeric(funcData$numSNPs)))
}

AllocateBinaryDosageWriteMemory43 <- function(funcData) {
  AllocateBinaryDosageWriteMemory41(funcData)
}

AllocateBinaryDosageWriteMemory44 <- function(funcData) {
  AllocateBinaryDosageWriteMemory42(funcData)
}

# Write binary dosage data at the end of the file
# Header has already been written
# funcData was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageFileData <- function(funcData, dosage, p0, p1, p2) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageData11, WriteBinaryDosageData12),
                    f2 <- c(WriteBinaryDosageData21, WriteBinaryDosageData22),
                    f3 <- c(WriteBinaryDosageData31, WriteBinaryDosageData32, WriteBinaryDosageData33, WriteBinaryDosageData34),
                    f4 <- c(WriteBinaryDosageData41, WriteBinaryDosageData42, WriteBinaryDosageData43, WriteBinaryDosageData44))
  return (writeFunc[[funcData$format]][[funcData$subformat]](funcData, dosage, p0, p1, p2))
}

WriteBinaryDosageData11 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 0))
}

WriteBinaryDosageData12 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(funcData$filename, p1, p2, funcData$usp1, funcData$usp2, 1))
}

WriteBinaryDosageData21 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData22 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(funcData$filename, p1, p2, funcData$usp1, funcData$usp2, 2))
}

WriteBinaryDosageData31 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData32 <- function(funcData, dosage, p0, p1, p2) {
  return (0)
}

WriteBinaryDosageData33 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData34 <- function(funcData, dosage, p0, p1, p2) {
  return (0)
}

WriteBinaryDosageData41 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData42 <- function(funcData, dosage, p0, p1, p2) {
  return (0)
}

WriteBinaryDosageData43 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageData(funcData$filename, dosage, funcData$usdosage, 2))
}

WriteBinaryDosageData44 <- function(funcData, dosage, p0, p1, p2) {
  return (0)
}

