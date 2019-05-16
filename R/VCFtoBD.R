#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param version The subversion of the format of the output
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
VCFtoBD <- function(vcfFile, bdFile, famFile = "", mapFile = "", format = 4L, version = 0L) {
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

  if (is.numeric(version) == FALSE && is.integer(version) == FALSE)
    stop("version must be an integer value")
  if (length(version) != 1)
    stop("version must be a single integer value")
  if (is.numeric(version) == TRUE) {
    if (floor(version) != version)
      stop("version must be an integer")
    version = floor(version)
  }
  if (version < 0 || version > 2)
    stop("version must be an integer value from 0 to 2")

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
  return (list());
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
GetVCFInfo <- function(vcfFile, index = TRUE) {
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
  vcfInfo$SNPs <- read.table(filename,
                             skip = vcfInfo$headersize,
                             colClasses = coltypes,
                             stringsAsFactors = FALSE)
  colnames(vcfInfo$SNPs) <- c("Chromosome", "Location", "SNPID",
                              "Reference", "Alternate", "Quality",
                              "Filter", "Info", "Format")
  return (vcfInfo)
}
