summarizevcfadditionalinfo <- function(x) {
  if (length(unique(x)) != 1)
    return (x)
  if (x[1] == '.')
    return (character(0))
  return (x[1])
}

readminimacinfofile <- function(filename) {
  addinfo <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  if (ncol(addinfo) != 13)
    stop("File dose not appear to be a minimac information file")
  if (all(colnames(addinfo) != c("SNP", "REF.0.", "ALT.1.", "ALT_Frq", "MAF",
                                 "AvgCall", "Rsq", "Genotyped", "LooRsq",
                                 "EmpR", "EmpRsq", "Dose0", "Dose1")))
    stop("File dose not appear to be a minimac information file")
  return(addinfo)
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
#' @param infofile Name of file name with information about the
#' imputation of the SNPs. This file is produced by minimac 3 and 4.
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
getvcfinfo <- function(filename,
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
  headersize <- -1L
  while (TRUE) {
    currentpos <- seek(con, origin = "current")
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
      begindata <- seek(con, origin = "current")
      break
    }
    headerlines <- headerlines + 1L
  }
  if (index == TRUE)
    headersize <- begindata
  close(con)

  if (all(x[1:9] == c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) == FALSE)
    stop("Column names incorrect")

  samples = data.frame(fid = rep("", length(x) - 9L),
                       sid = x[10:length(x)],
                       stringsAsFactors = FALSE)
  usesfid = FALSE

  coltypes = c("character", "integer", rep("character", 7),
               rep("NULL", nrow(samples)))
  snps <- read.table(filename,
                     skip = headerlines,
                     colClasses = coltypes,
                     stringsAsFactors = FALSE)
  colnames(snps) <- c("chromosome", "location", "snpid",
                      "reference", "alternate", "quality",
                      "filter", "info", "format")

  vcfinfo <- as.list(snps[,6:9])
  summaryinfo <- lapply(vcfinfo, summarizevcfadditionalinfo)
  datacolumns <- data.frame(numcolumns = rep(0L, length(summaryinfo$format)),
                            dosage = rep(0L, length(summaryinfo$format)),
                            genotypeprob = rep(0L, length(summaryinfo$format)),
                            genotype = rep(0L, length(summaryinfo$format)),
                            stringsAsFactors = FALSE)
  for (i in 1:length(summaryinfo$format)) {
    formatsplit <- unlist(strsplit(summaryinfo$format[i], split = ':'))
    datacolumns$numcolumns[i] <- length(formatsplit)
    datacolumns$dosage[i] <- match("DS", formatsplit)
    datacolumns$genotypeprob[i] <- match("GP", formatsplit)
    datacolumns$genotype[i] <- match("GT", formatsplit)
  }
  additionalinfo <- list(gzipped = gzipped,
                         headerlines = headerlines,
                         headersize = headersize,
                         quality = summaryinfo$quality,
                         filter = summaryinfo$filter,
                         info = summaryinfo$info,
                         format = summaryinfo$format,
                         datacolumns = datacolumns)
  class(additionalinfo) <- "vcf-info"
  rm(vcfinfo)
  rm(datacolumns)

  if (infofile == "") {
    snpinfo <- list()
  } else {
    minimacinfo <- readminimacinfofile(infofile)
    if (all(minimacinfo$SNP == snps$snpid) == TRUE &
        all(minimacinfo$REF.0. == snps$reference) == TRUE &
        all(minimacinfo$ALT.1. == snps$alternate) == TRUE) {
      snpinfo <- list(aaf = as.matrix(minimacinfo$ALT_Frq),
                      maf = as.matrix(minimacinfo$MAF),
                      avgcall = as.matrix(minimacinfo$AvgCall),
                      rsq = as.matrix(minimacinfo$Rsq))
    } else {
      stop("Imputation infromation file does not line up with VCF file")
    }
  }

  snps <- snps[,1:5]
  chr1 <- snps$chromosome[1]
  oneChr <- all(snps$chromosome == chr1)
  chrlocid <- paste(snps$chromosome, snps$location, sep = ":")
  vcfsnpformat1 <- all(snps$snpid == chrlocid)
  chrlocrefaltid <- paste(snps$chromosome, snps$location,
                          snps$reference, snps$alternate, sep = ":")
  vcfsnpformat2 <- all(snps$snpid == chrlocrefaltid)
  if (snpidformat == 0) {
    if (vcfsnpformat1 == TRUE) {
      snps$snpid <- chrlocid
      snpidformat <- 1L
    } else if (vcfsnpformat2 == TRUE) {
      snps$snpid <- chrlocrefaltid
      snpidformat <- 2L
    }
  } else if (snpidformat == 1) {
    if (vcfsnpformat2 == TRUE)
      stop ("snpidformat 1 specified but VCF file uses snpidformat 2")
    if (vcfsnpformat1 == FALSE)
      snps$snpid <- chrlocid
  } else if (snpidformat == 2) {
    if (vcfsnpformat2 == FALSE)
      snps$snpid <- chrlocrefaltid
  }

  if (index == TRUE) {
    datasize <- integer(nrow(snps))
    indices <- numeric(nrow(snps))
    if (gzipped == FALSE) {
      x <- GetLineLocations(filename)
      indices <- x[(headerlines + 1):(length(x) - 1)]
      for (i in 1:length(datasize))
        datasize[i] <- x[headerlines + i + 1] - x[headerlines + i]
    } else {
      con2 <- gzfile(filename, "r")
      currentpos <- seek(con2, 0)
      for (i in 1:headerlines)
        line <- readLines(con2, n = 1)
      currentpos <- 0
      for (i in 1:nrow(snps)) {
        indices[i] <- seek(con2)
        line <- readLines(con2, n = 1)
        currentpos <- seek(con2)
        datasize[i] <- currentPos - indices[i]
      }
      close(con2)
    }
  } else {
    datasize <- integer(0)
    indices <- numeric(0)
  }

  retval = list(filename = fqfilename,
                usesfid = usesfid,
                samples = samples,
                onechr = oneChr,
                snpidformat = snpidformat,
                snps = snps,
                snpinfo = snpinfo,
                datasize = datasize,
                indices = indices,
                additionalinfo = additionalinfo)
  class(retval) <- c("genetic-file-info", "vcf-file-info")
  return (retval)
}

#' Apply a function to all the SNPs in a VCF file
#'
#' Routine to apply a user provided function to all the
#' SNPs in a VCF file.
#'
#' @param vcfinfo Information about a VCF file returned from
#' getvcfinfo
#' @param func
#' The function to apply to each of the SNPs.
#'
#' @return
#' A list with the function result for each SNP
#'
#' @export
#'
#' @examples
#' # In work
vcfapply <- function(vcfinfo, func, ...) {
  if (is.na(match("vcf-file-info", class(vcfinfo))) == TRUE)
    stop("vcfinfo does not appear to contain information about a vcf file")

  retval <- vector("list", nrow(vcfinfo$snps))
  if (vcfinfo$additionalinfo$gzipped == FALSE)
    con <- file(vcfinfo$filename, "r")
  else
    con <- gzfile(vcfinfo$filename, "r")
  line <- readLines(con, n = vcfinfo$additionalinfo$headerlines)

  for (i in 1:nrow(vcfinfo$snps)) {
    line <- readLines(con, n = 1)
    x <- unlist(strsplit(line, "\t"))
    y <- unlist(strsplit(x[10:length(x)], ":"))
    if (length(vcfinfo$additionalinfo$datacolumns$dosage) == 1) {
      dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage
      gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob
      numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns
    } else {
      dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage[i]
      gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob[i]
      numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns[i]
    }
    if(is.na(dosagecol) == FALSE) {
      dosage <- as.numeric(y[seq(dosagecol, length(y) - numcolumns + dosagecol, numcolumns)])
    }
    if(is.na(gpcol) == FALSE) {
      gpstring <- y[seq(gpcol, length(y) - numcolumns + gpcol, numcolumns)]
      z <- unlist(strsplit(gpstring, ","))
      p0 <- as.numeric(z[seq(1, length(z) - 2, 3)])
      p1 <- as.numeric(z[seq(2, length(z) - 1, 3)])
      p2 <- as.numeric(z[seq(3, length(z), 3)])
    }
    retval[[i]] <- func(dosage = dosage,
                        p0 = p0,
                        p1 = p1,
                        p2 = p2
                        , ...)
  }
  close(con)
  return (retval)
}
