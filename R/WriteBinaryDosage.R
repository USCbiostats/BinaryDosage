#***************************************************************************#
#                                                                           #
#                       Writing Binary Dosage files                         #
#                                                                           #
#***************************************************************************#


#***************************************************************************#
#             Support functions for subject and SNP data                    #
#***************************************************************************#

# Formats 1, 2, and 3 all have a separate family and map file
# This routine saves the data frames in RDS format
WriteFamilyAndMapFiles <- function(filename, samples, snps) {
  saveRDS(samples, filename[2])
  saveRDS(snps, filename[3])
  return (md5 <- c(digest(samples, "md5"), digest(snps, "md5")))
}

# Find the groups value in the genetic file info. If it doesn't
# exist it returns the number of samples
FindGroups <- function(geneticfileinfo) {
  x <- match("groups", names(geneticfileinfo))
  if (is.na(x))
    return (geneticfileinfo$numSamples)
  return(geneticfileinfo$groups)
}

# Create the subject and family strings to write to the binary
# dosage header for format 4
SIDandFID4 <- function(funcData) {
  sid <- paste(funcData$samples$SID, collapse = '\t')
  if (funcData$usesFID == TRUE)
    fid <- paste(funcData$samples$FID, collapse = '\t')
  else
    fid <- ""
  return (list(sid = sid, fid = fid))
}

# Find a character vector in the SNP info data frame
# Returns the vector if found, otherwise it returns
# a character vector of length 0
# This is usually used to find aaf, maf, avgCall, and rsq
FindSNPInfoString <- function(toFind, snpInfo, numSNPs) {
  infocol <- match(toFind, colnames(snpInfo))
  if (is.na(infocol) == FALSE)
    return (snpInfo[,infocol])
  return ("")
}

# Find a numeric vector in the SNP info data frame or the
# bdoptions. If option is specified, a vector of zeros is
# return. If it is not found in the options, it
# returns the vector if found in the SNP info, otherwise
# it returns a numeric vector of length 0
FindSNPInfoNumeric <- function(toFind, snpInfo, numSNPs, bdoptions) {
  if (is.na(match(toFind, bdoptions)) == FALSE)
    return (numeric(numSNPs))
  infocol <- match(toFind, names(snpInfo))
  if (is.na(infocol) == FALSE)
    return (snpInfo[infocol])
  return (numeric(0))
}

# Find the SNP information needed for format 4
FindBDSNPInfo <- function(funcData, bdoptions) {
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
  return (list(snpid = snpid,
               chromosome = chr,
               location = loc,
               reference = ref,
               alternate = alt,
               aaf = aaf,
               maf = maf,
               avgCall = avgCall,
               rsq = rsq))
}

#***************************************************************************#
#                                                                           #
#                     Writing the Binary Dosage header                      #
#                                                                           #
#***************************************************************************#


# Writes the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
WriteBinaryDosageHeader <- function(format, subformat, filename, funcData, bdoptions) {
  writeHeaderFunc <- list(f1 <- c(WriteBinaryDosageHeader1, WriteBinaryDosageHeader1),
                          f2 <- c(WriteBinaryDosageHeader1, WriteBinaryDosageHeader1),
                          f3 <- c(WriteBinaryDosageHeader31, WriteBinaryDosageHeader31, WriteBinaryDosageHeader33, WriteBinaryDosageHeader34),
                          f4 <- c(WriteBinaryDosageHeader41, WriteBinaryDosageHeader41, WriteBinaryDosageHeader43, WriteBinaryDosageHeader44))
  WriteBinaryDosageBaseHeader(filename[1], format - 1, subformat - 1)
  return (writeHeaderFunc[[format]][[subformat]](filename, funcData, bdoptions))
}

WriteBinaryDosageHeader1 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  return (0)
}

WriteBinaryDosageHeader31 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3A(filename[1], funcData$numSamples)
  return (0)
}

WriteBinaryDosageHeader33 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2], 0)
  return (0)
}

WriteBinaryDosageHeader34 <- function(filename, funcData, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, funcData$samples, funcData$SNPs)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2], funcData$numSNPs)
  return (0)
}

WriteBinaryDosageHeader4 <- function(filename, funcData, bdoptions,
                                     headerEntries, offsets, numIndices) {
  subInfo <- SIDandFID4(funcData)
  snpInfo <- FindBDSNPInfo(funcData, bdoptions)
  WriteBinaryDosageHeader4A(filename[1],
                            headerEntries,
                            funcData$numSamples,
                            funcData$numSNPs,
                            FindGroups(funcData),
                            subInfo$sid[1],
                            subInfo$fid[1],
                            snpInfo$snpid[1],
                            snpInfo$chromosome[1],
                            snpInfo$location,
                            snpInfo$reference[1],
                            snpInfo$alternate[1],
                            snpInfo$aaf,
                            snpInfo$maf,
                            snpInfo$avgCall,
                            snpInfo$rsq,
                            offsets,
                            numIndices)
  return (0)
}

WriteBinaryDosageHeader41 <- function(filename, funcData, bdoptions) {
  headerEntries <- 8
  offsets <- c(seq(8L, 36L, 4L), 36L)
  return (WriteBinaryDosageHeader4(filename, funcData, bdoptions,
                                   headerEntries, offsets, 0))
}

WriteBinaryDosageHeader43 <- function(filename, funcData, bdoptions) {
  headerEntries <- 4
  offsets <- c(rep(-1L, 5), seq(8L, 20L, 4L))
  return (WriteBinaryDosageHeader4(filename, funcData, bdoptions,
                                   headerEntries, offsets, 0))
}

WriteBinaryDosageHeader44 <- function(filename, funcData, bdoptions) {
  headerEntries <- 4
  offsets <- c(rep(-1L, 5), seq(8L, 20L, 4L))
  return (WriteBinaryDosageHeader4(filename, funcData, bdoptions,
                                   headerEntries, offsets, funcData$numSNPs))
}

#***************************************************************************#
#                                                                           #
#                     Writing the Binary Dosage data                        #
#                                                                           #
#***************************************************************************#


# Allocates memory needed to write binary dosage files
# This is sufficient for all formats
AllocateBinaryDosageWriteMemory <- function(funcData) {
  return(list(filename = funcData$filename,
              format = funcData$format,
              subformat = funcData$subformat,
              headersize = funcData$headersize,
              snpnumber = 0L,
              datasize = integer(funcData$numSNPs),
              us = integer(2*funcData$numSamples)))
}

# Write binary dosage data at the end of the file
# Header has already been written
# funcData was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageData <- function(funcData, dosage, p0, p1, p2) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageData1, WriteBinaryDosageData2),
                    f2 <- c(WriteBinaryDosageData3, WriteBinaryDosageData4),
                    f3 <- c(WriteBinaryDosageData3, WriteBinaryDosageData5, WriteBinaryDosageData3, WriteBinaryDosageData6),
                    f4 <- c(WriteBinaryDosageData3, WriteBinaryDosageData5, WriteBinaryDosageData3, WriteBinaryDosageData6))
  return (writeFunc[[funcData$format]][[funcData$subformat]](funcData, dosage, p0, p1, p2))
}

WriteBinaryDosageData1 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageDataC(funcData$filename, dosage, funcData$us, 1))
}

WriteBinaryDosageData2 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(funcData$filename, p1, p2, funcData$us,  2))
}

WriteBinaryDosageData3 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryDosageDataC(funcData$filename, dosage, funcData$us, 3))
}

WriteBinaryDosageData4 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(funcData$filename, p1, p2, funcData$us, 3))
}

WriteBinaryDosageData5 <- function(funcData, dosage, p0, p1, p2) {
  snpnumber <- -1L
  return (WriteBinaryCompressed(funcData$filename,
                                dosage, p0, p1, p2,
                                snpnumber,
                                funcData$datasize,
                                funcData$us))
}

WriteBinaryDosageData6 <- function(funcData, dosage, p0, p1, p2) {
  return (WriteBinaryCompressed(funcData$filename,
                                dosage, p0, p1, p2,
                                funcData$snpnumber,
                                funcData$datasize,
                                funcData$us))
}

# Write binary dosage indices to  the file
# Header has already been written
# funcData was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageIndices <- function(funcData) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1),
                    f2 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1),
                    f3 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices2),
                    f4 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices2))
  return (writeFunc[[funcData$format]][[funcData$subformat]](funcData))
}

WriteBinaryDosageIndices1 <- function(funcData) {
  return (0)
}

WriteBinaryDosageIndices2 <- function(funcData) {
  return(WriteBinaryDosageIndicesC(funcData$filename, funcData$headersize, funcData$datasize))
}
