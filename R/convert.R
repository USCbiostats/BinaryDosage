#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
#' @importFrom prodlim row.match
#' @importFrom utils read.table
NULL

validatebdinput <- function(bdfiles,
                            format,
                            subformat,
                            snpidformat,
                            bdoptions) {
  if (is.character(bdfiles) == FALSE)
    stop("bdfiles must be a vector of characters")

  if (is.numeric(format) == FALSE && is.integer(format) == FALSE)
    stop("format must be an integer value")
  if (length(format) != 1)
    stop("format must be an integer vector of length 1")
  if (is.numeric(format) == TRUE) {
    if (floor(format) != format)
      stop("format must be an integer")
    format <- floor(format)
  }
  if (format < 1 || format > 4)
    stop("format must be an integer value from 1 to 4")

  if (is.numeric(subformat) == FALSE && is.integer(subformat) == FALSE)
    stop("subformat must be an integer value")
  if (length(subformat) != 1)
    stop("subformat must be an integer vector of length 1")
  if (is.numeric(subformat) == TRUE) {
    if (floor(subformat) != subformat)
      stop("subformat must be an integer")
    subformat <- floor(subformat)
  }
  if (subformat < 0 || subformat > 4)
    stop("subformat must be an integer value from 0 to 4")
  if (format < 3 && subformat > 2)
    stop("subformat must be an integer value from 0 to 2 for formats 1 and 2")

  if (format == 4 & length(bdfiles) != 1)
    stop("Only one output file name is needed when using format 4")
  if (format < 4 & length(bdfiles) != 3)
    stop("Three output file names are required when using formats 1, 2, and 3")
  if (is.na(match("", bdfiles)) == FALSE)
    stop("Output file names cannot be blank")

  if (is.numeric(snpidformat) == FALSE && is.integer(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be an integer vector of length 1")
  if (is.numeric(snpidformat) == TRUE) {
    if (floor(snpidformat) != snpidformat)
      stop("snpidformat must be an integer")
    snpidformat = floor(snpidformat)
  }
  if (snpidformat < -1 | snpidformat > 3)
    stop("snpidformat must be and integer from -1 to 3")

  if (is.character(bdoptions) == FALSE)
    stop("bdoptions must be a character array")
  if (length(bdoptions) > 0 & format != 4)
    stop("bdoptions can only be used with format 4")
  if (length(bdoptions) > 0) {
    if (any(is.na(match(bdoptions, c("aaf", "maf", "rsq")))) == TRUE)
      stop("Only valid bdoptions are aaf, maf, and rsq")
  }
  return(list(format = format,
              subformat = subformat))
}

###########################################################
#                  VCF to Binary Dosage                   #
###########################################################

#' Convert a VCF file to a binary dosage file
#'
#' Routine to read information from a VCF file and create
#' a binary dosage file. The function is designed to use
#' files return from the Michigan Imputation Server but will
#' run on other VCF files if they contain dosage and genetic
#' probabilities. Note: This routine can take a long time to
#' run if the VCF file is large.
#'
#' @param vcffiles A vector of file names.
#' The first is the name of the vcf file. The
#' second is name of the file that contains information
#' about the imputation of the SNPs. This file is produced
#' by minimac 3 and 4.
#' @param gz Indicator if VCF file is compressed using gzip.
#' Default value is FALSE.
#' @param bdfiles Vector of names of the output files.
#' The binary dosage file name is first. The family and
#' map files follow. For format 4, no family and map file
#' names are needed.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subformat of the format of the output
#' binary dosage file. A value of 1 or 3 indicates that only the
#' dosage value is saved. A value of 2 or 4 indicates
#' the dosage and genetic probabilities will be output. Values
#' of 3 or 4 are only allowed with formats 3 and 4. If a value
#' of zero if provided, and genetic probabilities are in the vcf
#' file, subformat 2 will be used for formats 1 and 2, and
#' subformat 4 will be used for formats 3 and 4. If the vcf file
#' does not contain genetic probabilities, subformat 1 will be
#' used for formats 1 and 2, and subformat 3 will be used for
#' formats 3 and 4. The default value is 0.
#' @param snpidformat The format that the SNP ID will be saved as.
#' -1 SNP ID not written
#' 0 - same as in the VCF file
#' 1 - chr:pos
#' 2 - chr:pos:ref:alt
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#' @param bdoptions Character array containing any of the following
#' value, "aaf", "maf", "rsq". The presence of any of these
#' values indicates that the specified values should be
#' calculates and stored in the binary dosage file. These values only
#' apply to format 4.
#'
#' @return
#' None
#' @export
#'
#' @examples
#' # Find the vcf file names
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
#' bdfiles <- tempfile()
#' # Convert the file
#' vcftobdlegacy(vcffiles = c(vcf1afile, vcf1ainfo), bdfiles = bdfiles)
#' # Verify the file was written correctly
#' bdinfo <- getbdinfo(bdfiles)
vcftobdlegacy <- function(vcffiles,
                    gz = FALSE,
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0,
                    bdoptions = character(0)) {
  message("vcftobdlegacy() uses legacy formats 1-4. ",
          "Consider using vcftobd() instead, which produces Format 5 files ",
          "with faster conversion and smaller file sizes.")
  if (missing(vcffiles) == TRUE)
    stop("No VCF file specified")

  if (missing(bdfiles) == TRUE)
    stop("No output files specified")

  validation <- validatebdinput(bdfiles = bdfiles,
                                format = format,
                                subformat = subformat,
                                snpidformat = snpidformat,
                                bdoptions = bdoptions)
  format <- validation$format
  subformat <- validation$subformat

  if (snpidformat == -1)
    readsnpformat = 0
  else
    readsnpformat = snpidformat
  vcfinfo <- getvcfinfo(vcffiles = vcffiles,
                        gz = gz,
                        index = FALSE,
                        snpidformat = readsnpformat)
  if (snpidformat == -1)
    vcfinfo$snpidformat = 1
  else
    vcfinfo$snpidformat = 0

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
  bdinfo <- getbdinfo(bdfiles = bdfiles)
  if (is.na(match("aaf", bdoptions)) == FALSE)
    updateaaf(bdinfo)
  if (is.na(match("maf", bdoptions)) == FALSE)
    updatemaf(bdinfo)
  if (is.na(match("rsq", bdoptions)) == FALSE)
    updatersq(bdinfo)
  ##return (0)
}

###########################################################
#          Binary Dosage Format 1 to Format 5            #
###########################################################

#' Convert a format 1 binary dosage file to Format 5
#'
#' Reads a format 1 binary dosage file (subformat 1 or 2) and converts it to
#' a Format 5 file pair. If the source file does not contain genotype
#' probabilities (subformat 1), those values are stored as missing in the
#' output.
#'
#' @param bdfiles Vector of three file names for the format 1 binary dosage
#'   file: the binary dosage file, the family file, and the map file.
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#'
#' @examples
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' bdfile  <- tempfile()
#' famfile <- tempfile()
#' mapfile <- tempfile()
#' suppressWarnings(
#'   vcftobdlegacy(vcffiles = vcf1afile,
#'                 bdfiles = c(bdfile, famfile, mapfile),
#'                 format = 1L)
#' )
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' bd1tobd5(bdfiles     = c(bdfile, famfile, mapfile),
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file)
bd1tobd5 <- function(bdfiles, bdose_file, bdinfo_file) {
  if (missing(bdfiles))
    stop("No binary dosage files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  old_bdinfo <- getbdinfo(bdfiles)

  fmt    <- old_bdinfo$additionalinfo$format
  subfmt <- old_bdinfo$additionalinfo$subformat

  if (fmt != 1L)
    stop("bdfiles must be a format 1 binary dosage file")
  if (!subfmt %in% c(1L, 2L))
    stop("bdfiles must be format 1 subformat 1 or 2")

  has_gp      <- subfmt == 2L
  n_samp      <- nrow(old_bdinfo$samples)
  n_snps      <- nrow(old_bdinfo$snps)
  indices     <- numeric(n_snps)
  current_pos <- 4.0  # bytes written after 4-byte magic header

  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  for (i in seq_len(n_snps)) {
    snp <- getsnp(old_bdinfo, i, dosageonly = FALSE)
    if (has_gp) {
      gp <- as.numeric(rbind(snp$p0, snp$p1, snp$p2))
    } else {
      gp <- rep(NA_real_, 3L * n_samp)
    }
    block      <- compress_snp_block(snp$dosage, gp)
    indices[i] <- current_pos
    writeBin(block, bdose_con)
    current_pos <- current_pos + length(block)
  }

  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = old_bdinfo$usesfid,
    samples        = old_bdinfo$samples,
    onechr         = old_bdinfo$onechr,
    snpidformat    = old_bdinfo$snpidformat,
    snps           = old_bdinfo$snps,
    snpinfo        = old_bdinfo$snpinfo,
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, bdinfo_file)

  invisible(NULL)
}

###########################################################
#          Binary Dosage Format 2 to Format 5            #
###########################################################

#' Convert a format 2 binary dosage file to Format 5
#'
#' Reads a format 2 binary dosage file (subformat 1 or 2) and converts it to
#' a Format 5 file pair. If the source file does not contain genotype
#' probabilities (subformat 1), those values are stored as missing in the
#' output.
#'
#' @param bdfiles Vector of three file names for the format 2 binary dosage
#'   file: the binary dosage file, the family file, and the map file.
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#'
#' @examples
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' bdfile  <- tempfile()
#' famfile <- tempfile()
#' mapfile <- tempfile()
#' suppressWarnings(
#'   vcftobdlegacy(vcffiles = vcf1afile,
#'                 bdfiles = c(bdfile, famfile, mapfile),
#'                 format = 2L)
#' )
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' bd2tobd5(bdfiles     = c(bdfile, famfile, mapfile),
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file)
bd2tobd5 <- function(bdfiles, bdose_file, bdinfo_file) {
  if (missing(bdfiles))
    stop("No binary dosage files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  old_bdinfo <- getbdinfo(bdfiles)

  fmt    <- old_bdinfo$additionalinfo$format
  subfmt <- old_bdinfo$additionalinfo$subformat

  if (fmt != 2L)
    stop("bdfiles must be a format 2 binary dosage file")
  if (!subfmt %in% c(1L, 2L))
    stop("bdfiles must be format 2 subformat 1 or 2")

  has_gp      <- subfmt == 2L
  n_samp      <- nrow(old_bdinfo$samples)
  n_snps      <- nrow(old_bdinfo$snps)
  indices     <- numeric(n_snps)
  current_pos <- 4.0

  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  for (i in seq_len(n_snps)) {
    snp <- getsnp(old_bdinfo, i, dosageonly = FALSE)
    if (has_gp) {
      gp <- as.numeric(rbind(snp$p0, snp$p1, snp$p2))
    } else {
      gp <- rep(NA_real_, 3L * n_samp)
    }
    block      <- compress_snp_block(snp$dosage, gp)
    indices[i] <- current_pos
    writeBin(block, bdose_con)
    current_pos <- current_pos + length(block)
  }

  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = old_bdinfo$usesfid,
    samples        = old_bdinfo$samples,
    onechr         = old_bdinfo$onechr,
    snpidformat    = old_bdinfo$snpidformat,
    snps           = old_bdinfo$snps,
    snpinfo        = old_bdinfo$snpinfo,
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, bdinfo_file)

  invisible(NULL)
}

###########################################################
#          Binary Dosage Format 3 to Format 5            #
###########################################################

#' Convert a format 3 binary dosage file to Format 5
#'
#' Reads a format 3 binary dosage file (subformat 1, 2, 3, or 4) and converts
#' it to a Format 5 file pair. If the source file does not contain genotype
#' probabilities (subformats 1 and 3), those values are stored as missing in
#' the output.
#'
#' @param bdfiles Vector of three file names for the format 3 binary dosage
#'   file: the binary dosage file, the family file, and the map file.
#'
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#'
#' @examples
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' bdfile  <- tempfile()
#' famfile <- tempfile()
#' mapfile <- tempfile()
#' suppressWarnings(
#'   vcftobdlegacy(vcffiles = vcf1afile,
#'                 bdfiles = c(bdfile, famfile, mapfile),
#'                 format = 3L)
#' )
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' bd3tobd5(bdfiles     = c(bdfile, famfile, mapfile),
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file)
bd3tobd5 <- function(bdfiles, bdose_file, bdinfo_file) {
  if (missing(bdfiles))
    stop("No binary dosage files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  old_bdinfo <- getbdinfo(bdfiles)

  fmt    <- old_bdinfo$additionalinfo$format
  subfmt <- old_bdinfo$additionalinfo$subformat

  if (fmt != 3L)
    stop("bdfiles must be a format 3 binary dosage file")
  if (!subfmt %in% c(1L, 2L, 3L, 4L))
    stop("bdfiles must be format 3 subformat 1, 2, 3, or 4")

  has_gp      <- subfmt %in% c(2L, 4L)
  n_samp      <- nrow(old_bdinfo$samples)
  n_snps      <- nrow(old_bdinfo$snps)
  indices     <- numeric(n_snps)
  current_pos <- 4.0

  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  for (i in seq_len(n_snps)) {
    snp <- getsnp(old_bdinfo, i, dosageonly = FALSE)
    if (has_gp) {
      gp <- as.numeric(rbind(snp$p0, snp$p1, snp$p2))
    } else {
      gp <- rep(NA_real_, 3L * n_samp)
    }
    block      <- compress_snp_block(snp$dosage, gp)
    indices[i] <- current_pos
    writeBin(block, bdose_con)
    current_pos <- current_pos + length(block)
  }

  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = old_bdinfo$usesfid,
    samples        = old_bdinfo$samples,
    onechr         = old_bdinfo$onechr,
    snpidformat    = old_bdinfo$snpidformat,
    snps           = old_bdinfo$snps,
    snpinfo        = old_bdinfo$snpinfo,
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, bdinfo_file)

  invisible(NULL)
}

###########################################################
#          Binary Dosage Format 4 to Format 5            #
###########################################################

#' Convert a format 4 binary dosage file to Format 5
#'
#' Reads a format 4 binary dosage file (subformat 1, 2, 3, or 4) and converts
#' it to a Format 5 file pair. If the source file does not contain genotype
#' probabilities (subformats 1 and 3), those values are stored as missing in
#' the output.
#'
#' @param bdfile File name of the format 4 binary dosage file. Format 4 stores
#'   all data in a single file.
#'
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#'
#' @return NULL (invisibly)
#'
#' @keywords internal
#'
#' @examples
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' bdfile <- tempfile()
#' suppressWarnings(
#'   vcftobdlegacy(vcffiles = vcf1afile,
#'                 bdfiles = bdfile,
#'                 format = 4L)
#' )
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' bd4tobd5(bdfile      = bdfile,
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file)
bd4tobd5 <- function(bdfile, bdose_file, bdinfo_file) {
  if (missing(bdfile))
    stop("No binary dosage file specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  old_bdinfo <- getbdinfo(bdfile)

  fmt    <- old_bdinfo$additionalinfo$format
  subfmt <- old_bdinfo$additionalinfo$subformat

  if (fmt != 4L)
    stop("bdfile must be a format 4 binary dosage file")
  if (!subfmt %in% c(1L, 2L, 3L, 4L))
    stop("bdfile must be format 4 subformat 1, 2, 3, or 4")

  has_gp      <- subfmt %in% c(2L, 4L)
  n_samp      <- nrow(old_bdinfo$samples)
  n_snps      <- nrow(old_bdinfo$snps)
  indices     <- numeric(n_snps)
  current_pos <- 4.0

  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  for (i in seq_len(n_snps)) {
    snp <- getsnp(old_bdinfo, i, dosageonly = FALSE)
    if (has_gp) {
      gp <- as.numeric(rbind(snp$p0, snp$p1, snp$p2))
    } else {
      gp <- rep(NA_real_, 3L * n_samp)
    }
    block      <- compress_snp_block(snp$dosage, gp)
    indices[i] <- current_pos
    writeBin(block, bdose_con)
    current_pos <- current_pos + length(block)
  }

  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = old_bdinfo$usesfid,
    samples        = old_bdinfo$samples,
    onechr         = old_bdinfo$onechr,
    snpidformat    = old_bdinfo$snpidformat,
    snps           = old_bdinfo$snps,
    snpinfo        = old_bdinfo$snpinfo,
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, bdinfo_file)

  invisible(NULL)
}

###########################################################
#        Update Binary Dosage Format 1-4 to Format 5     #
###########################################################

#' Update a binary dosage file to Format 5
#'
#' Reads a binary dosage file in format 1, 2, 3, or 4, detects the format
#' automatically, and converts it to a Format 5 file pair by calling the
#' appropriate conversion routine. If the source file does not contain genotype
#' probabilities, those values are stored as missing in the output.
#'
#' @param bdfiles Vector of file names for the binary dosage file. Format 4
#'   files require one file name. Formats 1, 2, and 3 require three file names:
#'   the binary dosage file, the family file, and the map file.
#'
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#'
#' @return NULL (invisibly)
#' @export
#'
#' @examples
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' bdfile <- tempfile()
#' suppressWarnings(
#'   vcftobdlegacy(vcffiles = vcf1afile,
#'                 bdfiles = bdfile,
#'                 format = 4L)
#' )
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' updatebd(bdfiles     = bdfile,
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file)
updatebd <- function(bdfiles, bdose_file, bdinfo_file) {
  if (missing(bdfiles))
    stop("No binary dosage files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  old_bdinfo <- getbdinfo(bdfiles)
  fmt        <- old_bdinfo$additionalinfo$format

  if (fmt == 1L)
    bd1tobd5(bdfiles = bdfiles, bdose_file = bdose_file, bdinfo_file = bdinfo_file)
  else if (fmt == 2L)
    bd2tobd5(bdfiles = bdfiles, bdose_file = bdose_file, bdinfo_file = bdinfo_file)
  else if (fmt == 3L)
    bd3tobd5(bdfiles = bdfiles, bdose_file = bdose_file, bdinfo_file = bdinfo_file)
  else if (fmt == 4L)
    bd4tobd5(bdfile = bdfiles[1L], bdose_file = bdose_file, bdinfo_file = bdinfo_file)
  else
    stop("bdfiles must be a format 1, 2, 3, or 4 binary dosage file")

  invisible(NULL)
}

###########################################################
#                  Gen to Binary Dosage                   #
###########################################################

#' Convert a gen file to a binary dosage file
#'
#' Routine to read information from a gen file and create
#' a binary dosage file. Note: This routine can take a long
#' time to run if the gen file is large.
#'
#' @param genfiles A vector of file names.
#' The first is the name of the gen file. The
#' second is name of the sample file that contains
#' the subject information.
#' @param snpcolumns Column numbers containing chromosome,
#' snpid, location, reference allele, alternate allele,
#' respectively. This must be an integer vector. All
#' values must be positive except for the chromosome.
#' The value for the chromosome may be -1 or -0.
#' -1 indicates that the chromosome value is passed to
#' the routine using the chromosome parameter.
#' 0 indicates that the chromosome value is in the snpid
#' and that the snpid has the format chromosome:other_data.
#' Default value is c(1L, 2L, 3L, 4L, 5L).
#' @param startcolumn Column number of first column with
#' genetic probabilities or dosages. Must
#' be an integer value. Default value is 6L.
#' @param impformat Number of genetic data values per
#' subject. 1 indicates dosage only, 2 indicates P(g=0)
#' and P(g=1) only, 3 indicates P(g=0), P(g=1), and
#' P(g=2). Default value is 3L.
#' @param chromosome Chromosome value to use if the
#' first value of the snpcolumns is equal to 0.
#' Default value is character().
#' @param header Indicators if the gen and sample files
#' have headers. If the gen file does not have a
#' header. A sample file must be included.
#' Default value is c(FALSE, TRUE).
#' @param gz Indicator if file is compressed using gzip.
#' Default value is FALSE.
#' @param sep Separator used in the gen file. Default
#' value is `"\t"`
#' @param bdfiles Vector of names of the output files.
#' The binary dosage file name is first. The family and
#' map files follow. For format 4, no family and map file
#' names are needed.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subformat of the format of the output
#' binary dosage file. A value of 1 or 3 indicates that only the
#' dosage value is saved. A value of 2 or 4 indicates
#' the dosage and genetic probabilities will be output. Values
#' of 3 or 4 are only allowed with formats 3 and 4. If a value
#' of zero if provided, and genetic probabilities are in the vcf
#' file, subformat 2 will be used for formats 1 and 2, and
#' subformat 4 will be used for formats 3 and 4. If the vcf file
#' does not contain genetic probabilities, subformat 1 will be
#' used for formats 1 and 2, and subformat 3 will be used for
#' formats 3 and 4. The default value is 0.
#' @param snpidformat The format that the SNP ID will be saved as.
#' -1 - SNP ID not written.
#' 0 - same as in the VCF file.
#' 1 - chr:pos.
#' 2 - chr:pos:ref:alt.
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#' @param bdoptions Character array containing any of the following
#' value, "aaf", "maf", "rsq". The presence of any of these
#' values indicates that the specified values should be
#' calculates and stored in the binary dosage file. These values only
#' apply to format 4.
#'
#' @return
#' None
#' @export
#'
#' @examples
#' # Find the gen file names
#' gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
#' gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
#' # Get temporary output file name
#' bdfiles <- tempfile()
#' # Convert the file
#' gentobd(genfiles = c(gen3afile, gen3asample),
#'         snpcolumns = c(0L, 2L:5L),
#'         bdfiles = bdfiles)
#' # Verify the file was written correctly
#' bdinfo <- getbdinfo(bdfiles = bdfiles)
gentobd <- function(genfiles,
                    snpcolumns = 1L:5L,
                    startcolumn = 6L,
                    impformat = 3L,
                    chromosome = character(),
                    header = c(FALSE, TRUE),
                    gz = FALSE,
                    sep = "\t",
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0L,
                    bdoptions = character(0)) {
  if (missing(genfiles) == TRUE)
    stop("No gen file specified")

  if (missing(bdfiles) == TRUE)
    stop("No output files specified")

  validation <- validatebdinput(bdfiles = bdfiles,
                                format = format,
                                subformat = subformat,
                                snpidformat = snpidformat,
                                bdoptions = bdoptions)
  format <- validation$format
  subformat <- validation$subformat

  if (snpidformat == -1)
    readsnpformat = 0L
  else
    readsnpformat = as.integer(snpidformat)
  geninfo <- getgeninfo(genfiles = genfiles,
                        snpcolumns = snpcolumns,
                        startcolumn = startcolumn,
                        impformat = impformat,
                        chromosome = chromosome,
                        header = header,
                        gz = gz,
                        index = FALSE,
                        snpidformat = readsnpformat,
                        sep = sep)
  if (snpidformat == -1)
    geninfo$snpidformat = 1
  else
    geninfo$snpidformat = 0

  if (subformat == 0) {
    if (geninfo$additionalinfo$format == 1L)
      subformat <- 1L
    else
      subformat <- 2L
  }
  WriteBinaryDosageHeader(format = format,
                          subformat = subformat,
                          filename = bdfiles,
                          genefileinfo = geninfo,
                          bdoptions = bdoptions)
  headerinfo <- ReadBinaryDosageHeader(filename = bdfiles)
  bdwriteinfo <- AllocateBinaryDosageWriteMemory(headerinfo = headerinfo)
  genapply(geninfo = geninfo,
           func = WriteBinaryDosageData,
           writeinfo = bdwriteinfo)
  WriteBinaryDosageIndices(writeinfo = bdwriteinfo)
  bdinfo <- getbdinfo(bdfiles = bdfiles)
  if (is.na(match("aaf", bdoptions)) == FALSE)
    updateaaf(bdinfo)
  if (is.na(match("maf", bdoptions)) == FALSE)
    updatemaf(bdinfo)
  if (is.na(match("rsq", bdoptions)) == FALSE)
    updatersq(bdinfo)
}

