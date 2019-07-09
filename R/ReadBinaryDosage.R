# Reads the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
ReadBinaryDosageHeader <- function(filename) {
  ReadHeaderFunc <- list(f1 <- c(ReadBinaryDosageHeader11, ReadBinaryDosageHeader12),
                         f2 <- c(ReadBinaryDosageHeader21, ReadBinaryDosageHeader22),
                         f3 <- c(ReadBinaryDosageHeader31, ReadBinaryDosageHeader32, ReadBinaryDosageHeader33, ReadBinaryDosageHeader34),
                         f4 <- c(ReadBinaryDosageHeader41, ReadBinaryDosageHeader42, ReadBinaryDosageHeader43, ReadBinaryDosageHeader44))
  bdformat <- ReadBinaryDosageBaseHeader(filename[1])
  return (ReadHeaderFunc[[bdformat$format]][[bdformat$subformat]](filename))
}

ReadFamAndMapFiles <- function(filename, format, subformat, headersize) {
  fqfilename <- normalizePath(filename[1], winslash = '/')
  subjects <- readRDS(filename[2])
  numsub <- nrow(subjects)
  if (all(subjects$FID == "") == TRUE)
    usesFID <- FALSE
  else
    usesFID <- TRUE
  snps <- readRDS(filename[3])
  numsnps <- nrow(snps)
  chr1 <- snps$Chromosome[1]
  onechr <- all(snps$Chromosome == chr1)
  chrLocID <- paste(snps$Chromosome, snps$Location, sep = ":")
  if(all(snps$SNPID == chrLocID) == TRUE)
    snpidformat <- 1
  else {
    chrLocRefAltID <- paste(snps$Chromosome, snps$Location,
                            snps$Reference, snps$Alternate, sep = ":")
    if (all(snps$SNPID == chrLocRefAltID) == TRUE)
      snpidformat <- 2
    else
      snpidformat <- 0
  }

  return (list(filename = fqfilename,
               format = format,
               subformat = subformat,
               headersize = headersize,
               numSamples = numsub,
               usesFID = usesFID,
               samples = subjects,
               onechr = onechr,
               snpidformat = snpidformat,
               numSNPs = numsnps,
               SNPs = snps))
}

ReadBinaryDosageHeader11 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 1, 1, 8))
}

ReadBinaryDosageHeader12 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 1, 2, 8))
}

ReadBinaryDosageHeader21 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 2, 1, 8))
}

ReadBinaryDosageHeader22 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 2, 2, 8))
}

ReadBinaryDosageHeader31 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 1, 12)
  additionalInfo = ReadBinaryDosageHeader3A(filename[1])
  if (additionalInfo$numsub != bdInfo$numSamples)
    stop("Subject file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader32 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 2, 12)
  additionalInfo = ReadBinaryDosageHeader3A(filename[1])
  if (additionalInfo$numsub != bdInfo$numSamples)
    stop("Subject file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader33 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 3, 72)
  additionalInfo = ReadBinaryDosageHeader3B(filename[1])
  if (digest(bdInfo$samples) != additionalInfo$md5[1])
    stop("Subject file does not line up with binary dosage file")
  if (digest(bdInfo$SNPs) != additionalInfo$md5[2])
    stop("Map file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader34 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 4, 72)
  additionalInfo = ReadBinaryDosageHeader3B(filename[1])
  if (digest(bdInfo$samples) != additionalInfo$md5[1])
    stop("Subject file does not line up with binary dosage file")
  if (digest(bdInfo$SNPs) != additionalInfo$md5[2])
    stop("Map file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader41 <- function(filename) {
  return (ReadBinaryDosageHeader4A(filename[1]))
}

ReadBinaryDosageHeader42 <- function(filename) {
  return (ReadBinaryDosageHeader4A(filename[1]))
}

ReadBinaryDosageHeader43 <- function(filename) {
  return (ReadBinaryDosageHeader4B(filename[1]))
}

ReadBinaryDosageHeader44 <- function(filename) {
  return (ReadBinaryDosageHeader4B(filename[1]))
}
