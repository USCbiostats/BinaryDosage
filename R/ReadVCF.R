SummarizeVCFAdditionalInfo <- function(x) {
  if (length(unique(x)) != 1)
    return (x)
  if (x[1] == '.')
    return (character(0))
  return (x[1])
}

#' Function to return information about a VCF file
#'
#' Function to return information about a VCF file.
#' This information is used by other routines to
#' allow for quicker extraction of values from the
#' file.
#'
#' @param filename Name of the VCF file.
#' @param gzipped Indicator if VCF file is compressed using gzip.
#' Default value is FALSE.
#' @param index Indicator if file should be indexed. This
#' allows for faster reading of the file. Indexing a gzipped
#' file does not speed up access. Default value is TRUE.
#' @param snpidformat The format that the SNP ID will be saved as.
#' 0 - same as in the VCF file
#' 1 - <chromosome>:<location>
#' 2 - <chromosome>:<location>:<reference allele>:<alternate allele>
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#'
#' @return List containing information about the VCF file
#' to include file name, subject IDs, and information about
#' the SNPs. Indices for faster reading will be included
#' if index is set to TRUE
#'
#' @export
#'
#' @examples
#' # Under construction
GetVCFInfo <- function(filename,
                       gzipped = FALSE,
                       index = TRUE,
                       snpidformat = 0L,
                       infofile = "") {
  if (missing(filename) == TRUE)
    stop("No VCF file specified")
  if (is.character(filename) == FALSE)
    stop("filename must be a character value")
  if (length(filename) != 1)
    stop("filename must be a single character value")
  if (filename == "")
    stop("No VCF file specified")

  if (is.logical(gzipped) == FALSE)
    stop("gzipped must be a logical value")
  if (length(gzipped) != 1)
    stop("gzipped must be a single logical value")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a single logical value")

  if (is.numeric(snpidformat) == FALSE)
    stop("snpidformat must be an intger value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be a singe integer value")
  if (floor(snpidformat) != snpidformat)
    stop("snpidformat must be an integer value")
  snpidformat <- as.integer(snpidformat)
  if (snpidformat < 0 || snpidformat > 2)
    stop("snpidformat must have a value of 0, 1, or 2")

  if (is.character(infofile) == FALSE)
    stop("infofile must be a chracter value")
  if (length(infofile) != 1)
    stop("infofile must be a single character value")

  if (gzipped == TRUE && index == TRUE)
    print("Indexing gzipped files is not recommended.")

  if (gzipped == FALSE) {
    con <- file(filename, "r")
  } else {
    con <- gzfile(filename, "r")
  }

  if (isOpen(con, "r") == FALSE)
    stop("Unable to open VCF file")
  fqfilename <- normalizePath(filename, winslash = '/')

  headerlines <- 1L
  while (TRUE) {
    currentPos <- seek(con, origin = "current")
    line <- readLines(con, n = 1)
    if (substr(line, 1, 1) != '#') {
      close(con)
      stop("Error process header")
    }
    if (substr(line, 2, 2) != '#') {
      x <- unlist(strsplit(line, "\t"))
      if (x[1] != "#CHROM") {
        close(con)
        stop("Error processing header")
      }
      x[1] = "CHROM"
      beginData <- seek(con, origin = "current")
      break
    }
    headerlines <- headerlines + 1L
  }
  close(con)

  if (all(x[1:9] == c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) == FALSE)
    stop("Column names incorrect")

  NumSamples = length(x) - 9L
  Samples = data.frame(FID = rep("", NumSamples),
                       SID = x[10:length(x)],
                       stringsAsFactors = FALSE)
  usesFID = FALSE

  coltypes = c("character", "integer", rep("character", 7),
               rep("NULL", NumSamples))
  SNPs <- read.table(filename,
                     skip = headerlines,
                     colClasses = coltypes,
                     stringsAsFactors = FALSE)
  colnames(SNPs) <- c("Chromosome", "Location", "SNPID",
                      "Reference", "Alternate", "quality",
                      "filter", "info", "format")
  VCFInfo <- as.list(SNPs[,6:9])
  SNPs <- SNPs[,1:5]

  numSNPs <- nrow(SNPs)
  chr1 <- SNPs$Chromosome[1]
  oneChr <- all(SNPs$Chromosome == chr1)
  chrLocID <- paste(SNPs$Chromosome, SNPs$Location, sep = ":")
  vcfSNPformat1 <- all(SNPs$SNPID == chrLocID)
  chrLocRefAltID <- paste(SNPs$Chromosome, SNPs$Location,
                          SNPs$Reference, SNPs$Alternate, sep = ":")
  vcfSNPformat2 <- all(SNPs$SNPID == chrLocRefAltID)
  if (snpidformat == 0) {
    if (vcfSNPformat1 == TRUE) {
      SNPs$SNPID <- chrLocID
      snpidformat <- 1L
    } else if (vcfSNPformat2 == TRUE) {
      SNPs$SNPID <- chrLocRefAltID
      snpidformat <- 2L
    }
  } else if (snpidformat == 1) {
    if (vcfSNPformat2 == TRUE)
      stop ("snpidformat 1 specified but VCF file uses snpidformat 2")
    if (vcfSNPformat1 == FALSE)
      SNPs$SNPID <- chrLocID
  } else if (snpidformat == 2) {
    if (vcfSNPformat2 == FALSE)
      SNPs$SNPID <- chrLocRefAltID
  }

  VCFInfo <- lapply(VCFInfo, SummarizeVCFAdditionalInfo)
  VCFInfo$dataColumns <- data.frame(numValues = rep(0L, length(VCFInfo$format)),
                                    dosage = rep(0L, length(VCFInfo$format)),
                                    genotypeProb = rep(0L, length(VCFInfo$format)),
                                    genotype = rep(0L, length(VCFInfo$format)),
                                    stringsAsFactors = FALSE)
  for (i in 1:length(VCFInfo$format)) {
    formatSplit <- unlist(strsplit(VCFInfo$format[i], split = ':'))
    VCFInfo$dataColumns$numValues[i] <- length(formatSplit)
    VCFInfo$dataColumns$dosage[i] <- match("DS", formatSplit)
    VCFInfo$dataColumns$genotypeProb[i] <- match("GP", formatSplit)
    VCFInfo$dataColumns$genotype[i] <- match("GT", formatSplit)
  }
  class(VCFInfo) <- "vcf-additional-info"

  if (index == TRUE) {
    datasize <- integer(numSNPs)
    Indices <- numeric(numSNPs)
    if (gzipped == FALSE) {
      x <- BinaryDosage:::GetLineLocations(filename)
      Indices <- x[(headerlines + 1):(length(x) - 1)]
      for (i in 1:length(datasize))
        datasize[i] <- x[headerlines + i + 1] - x[headerlines + i]
      headersize <- Indices[1]
    } else {
      con2 <- gzfile(filename, "r")
      currentPos <- seek(con2, 0)
      for (i in 1:headerlines)
        line <- readLines(con2, n = 1)
      headersize <- seek(con2)
      currentPos <- 0
      for (i in 1:numSNPs) {
        Indices[i] <- seek(con2)
        line <- readLines(con2, n = 1)
        currentPos <- seek(con2)
        datasize[i] <- currentPos - Indices[i]
      }
      close(con2)
    }
  } else {
    headersize <- -1L
    datasize <- integer(0)
    Indices <- numeric(0)
  }
  snpInfo <- list()

  retVal = list(filename = fqfilename,
                gzipped = gzipped,
                headerlines = headerlines,
                headersize = headersize,
                numSamples = NumSamples,
                usesFID = usesFID,
                samples = Samples,
                onechr = oneChr,
                snpidformat = snpidformat,
                numSNPs = numSNPs,
                SNPs = SNPs,
                snpInfo = snpInfo,
                datasize = datasize,
                indices = Indices,
                additionalInfo = VCFInfo)
  class(retVal) <- c("genetic-file-info", "vcf-file-info")
  return (retVal)
}

VCFApply <- function(vcfInfo, func, funcdata) {
  if (vcfInfo$gzipped == FALSE) {
    con <- file(vcfInfo$filename, "r")
    seek(con, sum(vcfInfo$indices[1]))
  } else {
    con <- gzfile(vcfInfo$filename, "r")
    line <- readLines(con, n = vcfInfo$headerlines)
  }

  for (i in 1:vcfInfo$numSNPs) {
    line <- readLines(con, n = 1)
    x <- unlist(strsplit(line, "\t"))
    y <- unlist(strsplit(x[10:length(x)], ":"))
    if (length(vcfInfo$additionalInfo$dataColumns$dosage) == 1) {
      dosageCol <- vcfInfo$additionalInfo$dataColumns$dosage
      gpCol <- vcfInfo$additionalInfo$dataColumns$genotypeProb
      numValues <- vcfInfo$additionalInfo$dataColumns$numValues
    } else {
      dosageCol <- vcfInfo$additionalInfo$dataColumns$dosage[i]
      gpCol <- vcfInfo$additionalInfo$dataColumns$genotypeProb[i]
      numValues <- vcfInfo$additionalInfo$dataColumns$numValues[i]
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
    func(funcdata, dosage, p0, p1, p2)
  }
  close(con)

  return (0)
}
