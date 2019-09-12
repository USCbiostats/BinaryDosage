gentobd <- function(genfile,
                    samplefile = "",
                    snpcolumns = 1L:5L,
                    startcolumn = 6L,
                    impformat = 3L,
                    chromosome = "",
                    header = FALSE,
                    gz = FALSE,
                    index = TRUE,
                    sep = '\t',
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0L,
                    bdoptions = character(0)) {
  if (missing(genfile) == TRUE)
    stop("No input file specified")
  if (is.character(genfile) == FALSE)
    stop("imp2file must be a character value")
  if (length(genfile) != 1)
    stop("imp2file must be a character vector of length 1")
  if (genfile == "")
    stop("No input file specified")

  if (is.character(samplefile) == FALSE)
    stop("samplefile must be a character value")
  if (length(samplefile) != 1)
    stop("samplefile must be a character vector of length 1")

  if (is.integer(snpcolumns) == FALSE)
    stop("snpcolumns must be an integer vector")
  if (length(snpcolumns) != 5)
    stop("snpcolumns must be an integer vector of length 5")
  if (min(snpcolumns[2:5]) < 1)
    stop("snpcolumns values other than chromosome must be positive integers")
  if (snpcolumns[1] < -1)
    stop("snpcolumns chromosome value must be -1, or a non-negatvie integer")

  if (is.integer(startcolumn) == FALSE)
    stop("startcolumn must be an integer value")
  if (length(startcolumn) != 1)
    stop("startcolumn must be an integer vector of length 1")
  if (startcolumn < 1)
    stop("startcolumn must be a positive integer")
  if (startcolumn <= max(snpcolumns))
    stop("startcolumn value must be larger than any value in snpcolumns")

  if (is.integer(impformat) == FALSE)
    stop("impformat must be an integer value")
  if (length(impformat) != 1)
    stop("impformat must be an integer vector of length 1")
  if (impformat < 1 | impformat > 3)
    stop("impformat must have a value of 1, 2, or 3")

  if (is.character(chromosome) == FALSE)
    stop("chromosome must be a character variable")
  if (length(chromosome) != 1)
    stop("chromosome must be a character vector of length 1")
  if (chromosome == "") {
    if (snpcolumns[1] == -1)
      stop("No chromosome column or value provided")
  } else {
    if (snpcolumns[1] > -1)
      stop("Both chromosome column and chromosome value provided")
  }

  if (is.logical(header) == FALSE)
    stop("header must be a logical value")
  if (length(header) != 1)
    stop("header must be a logical vector of length 1")
  if (header == FALSE) {
    if (samplefile == "")
      stop("File has no header and no sample file is provided")
  }

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a logical vector of length 1")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a logical vector of length 1")
  if (gz == TRUE & index == TRUE)
    stop("indexing of a gzipped file is not supported")

  if (is.integer(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be an integer vector of length 1")
  if (snpidformat < 0 | snpidformat > 3)
    stop("snpidformat must have a value of 0, 1, 2, or 3")

  if (is.character(sep) == FALSE)
    stop("sep must be a character value")
  if (length(sep) != 1)
    stop("sep must be a character vector of length 1")
  if (sep == "")
    stop("no value of sep provided")

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
    snpidformat <- as.integer(floor(snpidformat))
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

  if (is.character(bdoptions) == FALSE)
    stop("bdoptions must be a character array")
  if (length(bdoptions) > 0 & format != 4)
    stop("bdoptions can only be used with format 4")
  if (length(bdoptions) > 1) {
    if (any(is.na(match(bdoptions, c("aaf", "maf", "rsq")))) == TRUE)
      stop("Only valid bdoptions are aaf, maf, and rsq")
  }

  geninfo <- getgeninfo(genfiles = c(genfile, samplefile),
                        snpcolumns = snpcolumns,
                        startcolumn = startcolumn,
                        impformat = impformat,
                        chromosome = chromosome,
                        header = header,
                        gz = gz,
                        index = FALSE,
                        snpidformat = snpidformat,
                        sep = sep)

  if (subformat == 0) {
    if (anyNA(geninfo$additionalinfo$format) == 1)
      subformat <- 1
    else
      subformat <- 2
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
