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

# This routine sets up the indices when only the dosages are
# in the binary dosage file. This is simple because the size
# is fixed 2 bytes per subject per SNP
ReadIndices1 <- function(numSNPs, dosageoffset) {
  indices <- numeric(numSNPs)
  indices[1] <- dosageoffset
  for (i in 2:numSNPs)
    indices[i] = indices[i - 1] + 2 * numSNPs
  return (indices)
}

# This routine sets up the indices when reading formats 1 and 2
# when there is both dosage and genetic probability data. This
# is also simple because the size is fixed to 4 bytes per subject
# per SNP
ReadIndices2 <- function(numSNPs, dosageoffset) {
  indices <- numeric(numSNPs)
  indices[1] <- dosageoffset
  for (i in 2:numSNPs)
    indices[i] = indices[i - 1] + 4 * numSNPs
  return (indices)
}

# This routine sets up the indices when reading formats 3 and 4,
# subformat 2. The inidices are stored at the start of each SNP's
# dosage data. This is a time consuming operation.
ReadIndices3 <- function(numSNPs, dosageoffset) {
  indices <- numeric(numSNPs)
  indices[1] <- dosageoffset
  for (i in 2:numSNPs)
    indices[i] = indices[i - 1] + 4 * numSNPs
  return (indices)
}

# This routine sets up the indices when reading formats 3 and 4,
# subformat 4. The inidices are stored in the header are easy to
# read in.
ReadIndices4 <- function(filename, numSNPs, indexStart) {
  return (ReadBDIndicesS4(filename, numSNPs, indexStart))
}

Convert4HeaderToBDInfo <- function(filename, header, format, subformat) {
  SID <- unlist(strsplit(header$samples$sidstring, '\t'))
  if (header$samples$fidsize == 0) {
    usesFID <- FALSE
    FID = rep("", header$numsub)
  } else {
    usesFID <- TRUE
    FID <- unlist(strsplit(header$samples$fidstring, '\t'))
  }
  samples <- data.frame(FID, SID, stringsAsFactors = FALSE)

  if (header$snps$chrsize == 0) {
    onechr <- FALSE
    Chromosome <- rep("", header$numSNPs)
  } else {
    if (length(header$snps$chrstring) == 1) {
      onechr <- TRUE
      Chromosome <- rep(header$snps$chrstring, header$numSNPs)
    } else {
      onechr <- FALSE
      Chromosome <- unlist(strsplit(header$snps$chrstring, '\t'))
    }
  }
  if (length(header$snps$location) == 0)
    Location <- rep(0L, header$numSNPs)
  else
    Location <- header$snps$location
  if (header$snps$refsize == 0)
    Reference <- rep("", header$numSNPs)
  else
    Reference <- unlist(strsplit(header$snps$refstring, '\t'))
  if (header$snps$altsize == 0)
    Alternate <- rep("", header$numSNPs)
  else
    Alternate <- unlist(strsplit(header$snps$altstring, '\t'))
  if (header$snps$snpsize == 0)
    SNPID <- paste(Chromosome, Location, Reference, Alternate, sep = ':')
  else
    SNPID <- unlist(strsplit(header$snps$snpstring))
  SNPs <- data.frame(Chromosome, Location, SNPID, Reference, Alternate, stringsAsFactors = FALSE)

  snpInfoCol <- match(c("aaf", "maf", "avgcall", "rsq"), names(header$snps))
  snpInfoCol <- snpInfoCol[sapply(header$snps, function(x) length(x) != 0)[snpInfoCol]]
  snpinfo <- lapply(header$snps[snpInfoCol], matrix, nrow = header$numSNPs, ncol = header$numgroups)

  return (list(filename = normalizePath(filename[1], winslash = "/"),
               format = format,
               subformat = subformat,
               headersize = header$dosageoffset,
               numGroups = header$numgroups,
               groups = header$groups,
               numSamples = header$numsub,
               usesFID = usesFID,
               samples = samples,
               onechr = onechr,
               numSNPs = header$numSNPs,
               SNPs = SNPs,
               snpinfo = snpinfo))
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
  bdInfo$index <- ReadIndices4(filename[1], bdInfo$numSNPs, 72)
  bdInfo$headersize <- 72 + 4 * bdInfo$numSNPs
  return (bdInfo)
}

ReadBinaryDosageHeader41 <- function(filename) {
  header <- ReadBinaryDosageHeader4A(filename[1])
  return (Convert4HeaderToBDInfo(filename, header, 4, 1))
}

ReadBinaryDosageHeader42 <- function(filename) {
  header <- ReadBinaryDosageHeader4A(filename[1])
  return (Convert4HeaderToBDInfo(filename, header, 4, 2))
}

ReadBinaryDosageHeader43 <- function(filename) {
  header <- ReadBinaryDosageHeader4B(filename[1])
  bdInfo <- Convert4HeaderToBDInfo(filename, header, 4, 3)
  return (bdInfo)
}

ReadBinaryDosageHeader44 <- function(filename) {
  header <- ReadBinaryDosageHeader4B(filename[1])
  bdInfo <- Convert4HeaderToBDInfo(filename, header, 4, 4)
  bdInfo$index <- ReadIndices4(filename[1], bdInfo$numSNPs, header$indexoffset)
  return (bdInfo)
}

# Reads a SNP form the various formats
# of the binary dosage file.
ReadBinaryDosageData <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadHeaderFunc <- list(f1 <- c(ReadBinaryDosageData11, ReadBinaryDosageData12),
                         f2 <- c(ReadBinaryDosageData21, ReadBinaryDosageData22),
                         f3 <- c(ReadBinaryDosageData31, ReadBinaryDosageData32, ReadBinaryDosageData33, ReadBinaryDosageData34),
                         f4 <- c(ReadBinaryDosageData41, ReadBinaryDosageData42, ReadBinaryDosageData43, ReadBinaryDosageData44))
  return (ReadHeaderFunc[[bdInfo$format]][[bdInfo$subformat]](bdInfo, snp, d, p0, p1, p2, us))
}

ReadBinaryDosageData11 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 1)

}

ReadBinaryDosageData12 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataP1P2(bdInfo$filename, bdInfo$headersize, snp, d, p0, p1, p2, us, 2)
}

ReadBinaryDosageData21 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 3)

}

ReadBinaryDosageData22 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataP1P2(bdInfo$filename, bdInfo$headersize, snp, d, p0, p1, p2, us, 3)
}

ReadBinaryDosageData31 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 3)

}

ReadBinaryDosageData32 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  return (0)
}

ReadBinaryDosageData33 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 3)
}

ReadBinaryDosageData34 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  return (0)
}

ReadBinaryDosageData41 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 3)

}

ReadBinaryDosageData42 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  return (0)
}

ReadBinaryDosageData43 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  BinaryDosage:::ReadBinaryDosageDataC(bdInfo$filename, bdInfo$headersize, snp, d, us, 3)

}

ReadBinaryDosageData44 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  return (0)
}
