getimp2info <- function(imp2file,
                        samplefile = "",
                        snpcolumns = 1L:5L,
                        startcolumn = 6L,
                        impformat = 3L,
                        chromosome = "",
                        header = FALSE,
                        gzipped = FALSE,
                        index = TRUE,
                        snpidformat = 0L,
                        sep = '\t') {
  if (missing(imp2file) == TRUE)
    stop("No input file specified")
  if (is.character(imp2file) == FALSE)
    stop("imp2file must be a character value")
  if (length(imp2file) != 1)
    stop("imp2file must be a character vector of length 1")
  if (imp2file == "")
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

  if (is.logical(gzipped) == FALSE)
    stop("gzipped must be a logical value")
  if (length(gzipped) != 1)
    stop("gzipped must be a logical vector of length 1")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a logical vector of length 1")
  if (gzipped == TRUE & index == TRUE)
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


  if (header == TRUE) {
    if (gzipped == TRUE)
      filecon <- gzfile(imp2file)
    else
      filecon <- file(imp2file)
    headerline <- readLines(filecon, 1)
    close(filecon)
    headervalues <- unlist(strsplit(headerline, sep))
    if (length(headervalues) < startcolumn)
      stop("Number of values in header less than startcolumn")
    headervalues <- headervalues[startcolumn:length(headervalues)]
    if (length(headervalues) != 2 * floor(length(headervalues) / 2))
      stop("Odd number of values for family and subject ID")
    fid <- headervalues[seq(1,length(headervalues) - 1, 2)]
    iid <- headervalues[seq(2,length(headervalues), 2)]
    samples <- data.frame(fid = fid, sid = iid, stringsAsFactors = FALSE)
  } else {
    samples <- read.table(samplefile, header = TRUE, stringsAsFactors = FALSE)
    if (ncol(samples) == 1)
      stop("Error reading from sample file only one column found")
    samples <- samples[,1:2]
    samples[,1] <- as.character(samples[,1])
    samples[,2] <- as.character(samples[,2])
    colnames(samples) <- c("fid", "sid")
    if (samples[1,1] == "0" & samples[1,2] == "0")
      samples <- samples[2:nrow(samples),]
  }
  if (all(samples$fid == samples$sid)) {
    usesfid <- FALSE
    samples$fid <- ""
  } else {
    usesfid <- TRUE
  }

  coltypes = rep("NULL", impformat * nrow(samples) + (startcolumn - 1))
  if (snpcolumns[1] > 0)
    coltypes[snpcolumns[1]] <- "character"
  coltypes[snpcolumns[c(2, 4, 5)]] <- "character"
  coltypes[snpcolumns[3]] <- "integer"
  headersize <- 0
  if (header == TRUE)
    headersize <- 1
  snps <- read.table(imp2file,
                     skip = headersize,
                     colClasses = coltypes,
                     stringsAsFactors = FALSE)
  if (snpcolumns[1] == -1) {
    snps$chromosome <- chromosome
    snpcolumns[1] <- startcolumn
  } else if (snpcolumns[1] == 0) {
    chromosomeid <- unlist(strsplit(snps[,1], ':'))
    snpidentries <- length(chromosomeid) / nrow(snps)
    snps$chromosome <- chromosomeid[seq(1, length(chromosomeid) + 1 - snpidentries, snpidentries)]
    snpcolumns[1] <- startcolumn
  }
  x <- snpcolumns[2]
  snpcolumns[2] <- snpcolumns[3]
  snpcolumns[3] <- x
  snps <- snps[,order(order(snpcolumns))]
  colnames(snps) <- c("chromosome", "location", "snpid", "reference", "alternate")

  chr1 <- snps$chromosome[1]
  onechr <- FALSE
  if (all(snps$chromosome == chr1))
    onechr <- TRUE

  if (snpidformat == 0) {
    chrlocid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
    if (all(chrlocid == snps$snpid)) {
      snpidformat <- 1
    } else {
      chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location,
                                     snps$reference, snps$alternate, sep = ':'))
      if (all(chrlocrefaltid == snps$snpid)) {
        snpidformat <- 2
      } else {
        chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
        chrlocrefaltid <- unlist(paste(chrlocrefaltid, snps$reference,
                                       snps$alternate, sep = '_'))
        if (all(chrlocrefaltid == snps$snpid))
          snpidformat <- 3
      }
    }
  } else if (snpidformat == 1) {
    snps$snpid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
  } else if (snpidformat == 2) {
    snps$snpid <- unlist(paste(snps$chromosome, snps$location,
                               snps$reference, snps$alternate, sep = ':'))
  } else {
    chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
    snps$snpid <- unlist(paste(chrlocrefaltid, snps$reference,
                               snps$alternate, sep = '_'))
  }
  snpinfo <- list()
  if (index == FALSE) {
    datasize <- integer(0)
    indices <- numeric(0)
  } else {
    datasize <- integer(nrow(snps))
    indices <- numeric(nrow(snps))
    if (gzipped == FALSE) {
      x <- GetLineLocations(imp2file)
      if (header == TRUE)
        headerlines <- 1
      else
        headerlines <- 0
      indices <- x[(headerlines + 1):(length(x) - 1)]
      for (i in 1:length(datasize))
        datasize[i] <- x[headerlines + i + 1] - x[headerlines + i]
    } else {
      con2 <- gzfile(imp2file, "r")
      currentpos <- seek(con2, 0)
      if (header == TRUE)
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
  }

  additionalinfo <- list(gzipped = gzipped,
                         headersize = headersize,
                         format = impformat,
                         startcolumn = startcolumn)
  class(additionalinfo) <- "gen-info"

  retval <- list(filename = normalizePath(imp2file, winslash = '/'),
                 usesfid = usesfid,
                 samples = samples,
                 onechr = onechr,
                 snpidformat = snpidformat,
                 snps = snps,
                 snpinfo = snpinfo,
                 datasize = datasize,
                 indices = indices,
                 additionalinfo = additionalinfo)
  class(retval) <- "genetic-info"

  return(retval)
}

genapply <- function(geninfo, func, ...) {
  if (is.na(match("genetic-info", class(geninfo))) == TRUE)
    stop("geninfo is not of class genetic-info")
  if (is.na(match("vcf-info", class(vcfinfo$additionalinfo))) == TRUE)
    stop("geninfo does not appear to contain information about a gen file")

  retval <- vector("list", nrow(geninfo$snps))
  if (geninfo$additionalinfo$gzipped == FALSE)
    con <- file(geninfo$filename, "r")
  else
    con <- gzfile(geninfo$filename, "r")
  line <- readLines(con, n = geninfo$additionalinfo$headerlines)

  dosage <- numeric(nrow(geninfo$samples))
  p0 <- numeric(nrow(geninfo$samples))
  p1 <- numeric(nrow(geninfo$samples))
  p2 <- numeric(nrow(geninfo$samples))
  for (i in 1:nrow(vcfinfo$snps)) {
    line <- readLines(con, n = 1)
    x <- unlist(strsplit(line, "\t"))
    y <- unlist(strsplit(x[geninfo$additionalinfo$startcolumn:length(x)], ":"))
    if (geninfo$additionalinfo$format == 1) {
      dosage[1:length(dosage)] <- as.character(y)
    } else if (geninfo$additionalinfo$format == 2) {
      p0[1:length(p0)] <- as.character(y[seq(1, length(y) - 1, 2)])
      p1[1:length(p1)] <- as.character(y[seq(2, length(y), 2)])
      p2[1:length(p2)] <- ifelse(p0 + p1 < 1, 1 - p0 - p1, 0.)
      dosage[1:length(dosage)] <- ifelse(p1 + p2 + p2 < 2, p1 + p2 + p2, 2.)
    } else {
      p0[1:length(p0)] <- as.character(y[seq(1, length(y) - 2, 3)])
      p1[1:length(p1)] <- as.character(y[seq(2, length(y) - 1, 3)])
      p2[1:length(p2)] <- as.character(y[seq(3, length(y), 3)])
      dosage[1:length(dosage)] <- ifelse(p1 + p2 + p2 < 2, p1 + p2 + p2, 2.)
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
