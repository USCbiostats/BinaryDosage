#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
NULL

#' Function to return information about a VCF file
#'
#' Function to return information about a VCF file.
#' This information is used by other routines to
#' allow for quicker extraction of values from the
#' file.
#'
#' @param filename Name of the VCF file.
#' @param gz Indicator if VCF file is compressed using gzip.
#' Default value is FALSE.
#' @param index Indicator if file should be indexed. This
#' allows for faster reading of the file. WARNING: this is
#' not guaranteed to work on systems running Windows. Default
#' value is TRUE.
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
GetVCFInfo <- function(filename, gz = FALSE, index = TRUE) {
  if (missing(filename) == TRUE)
    stop("No VCF file specified")
  if (is.character(filename) == FALSE)
    stop("filename must be a character value")
  if (length(filename) != 1)
    stop("filename must be a single character value")
  if (filename == "")
    stop("No VCF file specified")

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a single logical value")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a single logical value")

  headersize <- 1L
  if (gz == TRUE && .Platform$OS.type == "windows" && index == TRUE)
    print("Indexing gz files in R on Windows is not recommended.")
  if (gz == FALSE) {
    con <- file(filename, "r")
  } else {
    con <- gzfile(filename, "r")
  }

  while (TRUE) {
    currentPos <- seek(con, origin = "current")
    line <- readLines(con, n = 1)
    if (substr(line, 1, 1) != '#') {
      close(con)
      return (-1L)
    }
    if (substr(line, 2, 2) != '#') {
      if (substr(line, 2, 6) != "CHROM") {
        close(con)
        return (-1L)
      }
      x <- unlist(strsplit(line, "\t"))
      beginData <- seek(con, origin = "current")
      break
    }
    headersize <- headersize + 1L
  }
  close(con)

  NumSamples = length(x) - 9L
  coltypes = c("character", "integer", rep("character", 7),
               rep("NULL", NumSamples))
  SNPs <- read.table(filename,
                     skip = headersize,
                     colClasses = coltypes,
                     stringsAsFactors = FALSE)
  colnames(SNPs) <- c("Chromosome", "Location", "SNPID",
                      "Reference", "Alternate", "Quality",
                      "Filter", "Info", "Format")
  numSNPs <- nrow(SNPs)
  chr1 <- SNPs$Chromosome[1]
  oneChr <- all(SNPs$Chromosome == chr1)
  chrLocID <- paste(SNPs$Chromosome, SNPs$Location, sep = ":")
  formattedID <- all(SNPs$SNPID == chrLocID)
  noQuality <- all(SNPs$Quality == '.')
  noInfo <- all(SNPs$Info == '.')
  format1 <- SNPs$Format[1]
  oneFormat <- all(SNPs$Format == format1)
  if (oneFormat) {
    formatSplit <- unlist(strsplit(format1, split = ':'))
    dosageColumn <- rep(match("DS", formatSplit), numSNPs)
    geneProbColumn <- rep(match("GP", formatSplit), numSNPs)
    genotypeColumn <- rep(match("GT", formatSplit), numSNPs)
  } else {
    dosageColumn <- integer(numSNPs)
    geneProbColumn <- integer(numSNPs)
    genotypeColumn <- integer(nuumSNPs)
    for (i in 1:numSNPs) {
      formatSplit <- unlist(strsplit(SNPs$Format[i], split = ':'))
      dosageColumn[i] <- match("DS", formatSplit)
      geneProbColumn[i] <- match("GP", formatSplit)
      genotypeColumn[i] <- match("GT", formatSplit)
    }
  }
  valueColumns <- data.frame(dosage = dosageColumn, geneProb = geneProbColumn, genotype = genotypeColumn)

  Indices <- numeric(numSNPs)
  if (gz == FALSE) {
    con2 <- file(filename, "r")
  } else {
    con2 <- gzfile(filename, "r")
  }

  currentPos <- seek(con2, 0)
  for (i in 1:headersize)
    line <- readLines(con2, n = 1)
  currentPos <- 0
  for (i in 1:numSNPs) {
    Indices[i] = seek(con2) - currentPos
    currentPos = seek(con2)
    line <- readLines(con2, n = 1)
  }
  close(con2)

  return (list(filename = filename,
               headersize = headersize,
               NumSamples = NumSamples,
               usesFID = FALSE,
               Samples = data.frame(FID = rep("", NumSamples),
                                    SID = x[10:length(x)],
                                    stringsAsFactors = FALSE),
               Indices = Indices,
               onechr = oneChr,
               formattedID = formattedID,
               noQuality = noQuality,
               noInfo = noInfo,
               oneFormat = oneFormat,
               NumSNPs = numSNPs,
               valueColumns = valueColumns,
               SNPs = SNPs))
}

#' Function to convert a VCF file to a binary dosage file
#'
#' Function to read information from a VCF file and create
#' a binary dosage file. The function is designed to use
#' files return from the Michigan Imputation Server but will
#' run on other VCF files if they contain dosage and genetic
#' probabilities.
#'
#' @param vcfFile Name of VCF file
#' @param bdFile Name of the output binary dosage file
#' @param famFile Name of output family information file.
#' This is only used if the output format is 1, 2, or 3.
#' The default value is ""
#' @param mapFile Name of the output map information file.
#' This is only used if the output format is 1, 2, or 3.
#' The default value is ""
#' @param gz Indicatator if vcf file is compressed using gzip.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subsubformat of the format of the output
#' binary dosage file. A value of 1 indicates that only the
#' dosage value is saved. A value greater than 1 indicates
#' the dosage and genetic probabilities will be output. A
#' value of 0 indicates the genetic probabilities will be
#' output if they are contained in the VCF file. Values of
#' 1 and 2 are currently supported for all formats.
#' The default value is 0.
#'
#' @return
#' A list containing information about the binary dosage file.
#' This is the same list returned from GetBDoseInfo. See
#' GetBDoseInfo for more information.
#' @export
#'
#' @examples
#' # Under construnction
VCFtoBD <- function(vcfFile, bdFile, famFile = "",
                    mapFile = "", gz = FALSE,
                    format = 4L, subformat = 0L) {
  if (missing(vcfFile) == TRUE)
    stop("No VCF file specified")
  if (is.character(vcfFile) == FALSE)
    stop("vcfFile must be a character value")
  if (length(vcfFile) != 1)
    stop("vcfFile must be a single character value")
  if (vcfFile == "")
    stop("No VCF file specified")

  if (missing(bdFile) == TRUE)
    stop("No output file specified")
  if (is.character(bdFile) == FALSE)
    stop("bdFile must be a character value")
  if (length(bdFile) != 1)
    stop("bdFile must be a single character value")
  if (bdFile == "")
    stop("No output file specified")

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a single logical value")

  if (is.numeric(format) == FALSE && is.integer(format) == FALSE)
    stop("format must be an integer value")
  if (length(format) != 1)
    stop("format must be a single integer value")
  if (is.numeric(format) == TRUE) {
    if (floor(format) != format)
      stop("format must be an integer")
    format = floor(format)
  }
  if (format < 1 || format > 4)
    stop("format must be an integer value from 1 to 4")

  if (is.numeric(subformat) == FALSE && is.integer(subformat) == FALSE)
    stop("subformat must be an integer value")
  if (length(subformat) != 1)
    stop("subformat must be a single integer value")
  if (is.numeric(subformat) == TRUE) {
    if (floor(subformat) != subformat)
      stop("subformat must be an integer")
    subformat = floor(subformat)
  }
  if (subformat < 0 || subformat > 4)
    stop("subformat must be an integer value from 0 to 4")
  if (format < 3 && subformat > 2)
    stop("subformat must be an integer value from 0 to 2 for formats 1 and 2")

  if (format < 4) {
    if (is.character(famFile) == FALSE)
      stop("famFile must be a character value")
    if (length(famFile) != 1)
      stop("famFile must be a single character value")
    if (famFile == "") {
      errormessage <- paste("famFile must be specified for format", format)
      stop(errormessage)
    }
    if (is.character(mapFile) == FALSE)
      stop("mapFile must be a character value")
    if (length(mapFile) != 1)
      stop("mapFile must be a single character value")
    if (mapFile == "") {
      errormessage <- paste("mapFile must be specified for format", format)
      stop(errormessage)
    }
  } else {
    if (is.character(famFile) == FALSE ||
        length(famFile) != 1 || famFile != "")
      stop("Value or famFile specified for format 4")
    if (is.character(mapFile) == FALSE ||
        length(mapFile) != 1 || mapFile != "")
      stop("Value or mapFile specified for format 4")
  }

  vcfInfo <- GetVCFInfo(vcfFile, gz = gz, index = FALSE)
  WriteBinaryDosageHeader(bdFile, format, subformat)
  WriteBDFamilyFile(bdFile, famFile, vcfInfo$Samples, format, subformat)
  return (vcfInfo)
}

#' Function to read information about a VCF file
#'
#' Function reads information about a VCf file. This
#' information is used to efficiently read dosages and
#' genetic probabilities from the VCF file using the
#' GetSNPValues function
#'
#' @param vcfFile Name of the VCF file
#' @param index Indicator if the file should be indexed
#' for quicker reader. This is not needed if only the
#' list of subjects and SNPs is needed. The default value
#' is TRUE.
#'
#' @return List of subjects and SNPs in VCF file. Indices
#' for reading the file are included if value index = TRUE
#' was passed to function.
#' @export
#'
#' @examples
#' # Under construction
GetVCFInfoC <- function(vcfFile, index = TRUE) {
  if (missing(vcfFile) == TRUE)
    stop("No VCF file specified")
  if (is.character(vcfFile) == FALSE)
    stop("vcfFile must be a character value")
  if (length(vcfFile) != 1)
    stop("vcfFile must be a single character value")
  if (vcfFile == "")
    stop("No VCF file specified")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a single logical value")

  vcfInfo <- OpenVCFFile(vcfFile)
  if ("genetic-file-info" %in% class(vcfInfo) == FALSE)
    stop("Error reading VCF information")
  coltypes <- c("character", "integer", rep("character", 7),
                rep("NULL", vcfInfo$NumSamples))
  vcfInfo$SNPs <- read.table(vcfFile,
                             skip = vcfInfo$headersize,
                             colClasses = coltypes,
                             stringsAsFactors = FALSE)
  colnames(vcfInfo$SNPs) <- c("Chromosome", "Location", "SNPID",
                              "Reference", "Alternate", "Quality",
                              "Filter", "Info", "Format")
  return (vcfInfo)
}

WriteBDFamilyFile <- function(bdFile, famFile, famInfo, format, subformat) {
  if (format > 3)
    return (1)
  if (all(famInfo[,"FID"] == "") == TRUE)
    usesFamilyID = FALSE
  else
    usesFamilyID = TRUE
  numSub = nrow(famInfo)
  famList <- list(numsub = numSub,
                  usesfamid = usesFamilyID,
                  md5 = digest::digest(famInfo, "md5"),
                  subjects = famInfo)
  class(famList) <- "famfile"
  saveRDS(famList, famFile)
  return (0)
}

WriteBDMapFile <- function(bdFile, mapFile, snpInfo, format, subformat) {
  if (format > 3)
    return (1)

  numSNPs <- nrow(snpInfo)
  chr1 <- snpInfo$Chromosome[1]
  oneChr <- all(snpInfo$Chromosome == chr1)
  chrLocID <- paste(snpInfo$Chromosome, snpInfo$Location, sep = ":")
  formattedID <- all(snpInfo$SNPID == chrLocID)
  noQuality <- all(snpInfo$Quality == '.')
  noInfo <- all(snpInfo$Info == '.')
  format1 <- snpInfo$Format[1]
  oneFormat <- all(snpInfo$Format == format1)
  md5 <- digest::digest(snpInfo)

  mapList <- list(numsnps = numSNPs,
                  onechr = oneChr,
                  formattedID = formattedID,
                  noQuality = noQuality,
                  noInfo = noInfo,
                  oneFormat = oneFormat,
                  md5 = md5,
                  snps = snpInfo)
  class(mapList) <- "mapfile"
  saveRDS(mapList, mapFile)
  return (0)
}
