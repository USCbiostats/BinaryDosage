#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
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
