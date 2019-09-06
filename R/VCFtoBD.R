#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
#' @importFrom prodlim row.match
NULL

#' Function to convert a VCF file to a binary dosage file
#'
#' Function to read information from a VCF file and create
#' a binary dosage file. The function is designed to use
#' files return from the Michigan Imputation Server but will
#' run on other VCF files if they contain dosage and genetic
#' probabilities.
#'
#' @param vcffile Name of VCF file
#' @param vcfinfofile (Optional) Name of information file associated
#' with the vcf file. Default value "".
#' @param gz Indicator if vcf file in compressed using gzip.
#' The default value is FALSE.
#' @param bdfiles Vector of names of the output files.
#' The binary dosage file name is first. The family and
#' map files follow. For format 4, no family and map file
#' names are needed.
#' @param snpidformat Format to store the snp id in. Format 1
#' is chromosome:location. Format 2 is
#' chromosome:location:reference:alternate. Format 0 indicates
#' to use the format in the vcf file. 0 is the default.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subsubformat of the format of the output
#' binary dosage file. A value of 1 or 3 indicates that only the
#' dosage value is saved. A value of 2 or 4 indicates
#' the dosage and genetic probabilities will be output. Values
#' of 3 or 4 are only allowed with formats 3 and 4. If a value
#' of zero if provided, and genetic probabilites are in the vcf
#' file, subformat 2 will be used for formats 1 and 2, and
#' subformat 4 will be used for formats 3 and 4. If the vcf file
#' does not contain genetic probabilities, subformat 1 will be
#' used for formats 1 and 2, and subformat 3 will be used for
#' formats 3 and 4. The default value is 0.
#' @param bdoptions Character array containg any of the following
#' value, "aaf", "maf", "rsq". The presence of any of these
#' values indicates that the specified values should be
#' calculates and stored in the bdosage file. These values only
#' apply to format 4.
#'
#' @return
#' A list containing information about the binary dosage file.
#' This is the same list returned from BDInfo. See
#' GetBDoseInfo for more information.
#' @export
#'
#' @examples
#' # Under construnction
vcftobd <- function(vcffile,
                    vcfinfofile = "",
                    gz = FALSE,
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0,
                    bdoptions = character(0)) {
  if (missing(vcffile) == TRUE)
    stop("No VCF file specified")
  if (is.character(vcffile) == FALSE)
    stop("vcffile must be a character value")
  if (length(vcffile) != 1)
    stop("vcffile must be a character array of length 1")
  if (vcffile == "")
    stop("No VCF file specified")

  if (is.character(vcfinfofile) == FALSE)
    stop("vcfinfofile must be a character value")
  if (length(vcfinfofile) != 1)
    stop("vcfinfofile must be a character array of length 1")

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

  if (is.numeric(snpidformat) == FALSE && is.integer(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be a single integer value")
  if (is.numeric(snpidformat) == TRUE) {
    if (floor(snpidformat) != snpidformat)
      stop("snpidformat must be an integer")
    snpidformat = floor(snpidformat)
  }
  if (snpidformat < 0 || snpidformat > 2)
    stop("snpidformat must be an integer value from 0 to 2")

  if (missing(bdfiles) == TRUE)
    stop("No output files specified")
  if (is.character(bdfiles) == FALSE)
    stop("Output file names must be a character values")
  if (format == 4 & length(bdfiles) != 1)
    stop("Only one file name is needed when using format 4")
  if (format < 4 & length(bdfiles) != 3)
    stop("Three file names are required when using formats 1, 2, and 3")
  if (is.na(match("", bdfiles)) == FALSE)
    stop("Output file names cannot be blank")

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a single logical value")

  if (is.character(bdoptions) == FALSE)
    stop("bdoptions must be a character array")
  if (length(bdoptions) > 0 & format != 4)
    stop("bdoptions can only be used with format 4")
  if (length(bdoptions) > 1) {
    if (any(is.na(match(bdoptions, c("aaf", "maf", "rsq")))) == TRUE)
      stop("Only valid bdoptions are aaf, maf, and rsq")
  }

  vcfinfo <- getvcfinfo(filename = vcffile,
                        gz = gz,
                        index = FALSE,
                        snpidformat = snpidformat)
  if (subformat == 0) {
    if (anyNA(vcfinfo$additionalinfo$datacolumns$genotypeprob) == TRUE)
      subformat <- 1
    else
      subformat <- 2
  }
  WriteBinaryDosageHeader(format = format,
                          subformat = subformat,
                          filename = bdfiles,
                          genefileinfo = vcfinfo,
                          bdoptions = bdoptions)
  headerinfo <- ReadBinaryDosageHeader(filename = bdfiles)
  bdwriteinfo <- AllocateBinaryDosageWriteMemory(headerinfo = headerinfo)
  vcfapply(vcfinfo = vcfinfo,
           func = WriteBinaryDosageData,
           writeinfo = bdwriteinfo)
  WriteBinaryDosageIndices(writeinfo = bdwriteinfo)
  bdinfo <- getbdinfo(bdfilenames = bdfiles)
  if (is.na(match("aaf", bdoptions)) == FALSE)
    updateaaf(bdinfo)
  if (is.na(match("maf", bdoptions)) == FALSE)
    updatemaf(bdinfo)
  if (is.na(match("rsq", bdoptions)) == FALSE)
    updatersq(bdinfo)
}
