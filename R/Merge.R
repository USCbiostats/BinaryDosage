mergesnpinfo <- function (mergedinfo,
                          numgroups,
                          geneticinfo,
                          snpsbtom,
                          snpsmtob) {
  hasaaf <- FALSE
  hasmaf <- FALSE
  hasavgcall <- FALSE
  hasrsq <- FALSE

  if (is.na(match("aaf", names(geneticinfo[[1]]$snpinfo))) == FALSE)
    hasaaf <- TRUE
  if (is.na(match("maf", names(geneticinfo[[1]]$snpinfo))) == FALSE)
    hasmaf <- TRUE
  if (is.na(match("avgcall", names(geneticinfo[[1]]$snpinfo))) == FALSE)
    hasavgcall <- TRUE
  if (is.na(match("rsq", names(geneticinfo[[1]]$snpinfo))) == FALSE)
    hasrsq <- TRUE
  for (i in 2:length(geneticinfo)) {
    if (is.na(match("aaf", names(geneticinfo[[i]]$snpinfo))) == FALSE)
      hasaaf <- TRUE
    if (is.na(match("maf", names(geneticinfo[[i]]$snpinfo))) == FALSE)
      hasmaf <- TRUE
    if (is.na(match("avgcall", names(geneticinfo[[i]]$snpinfo))) == FALSE)
      hasavgcall <- TRUE
    if (is.na(match("rsq", names(geneticinfo[[i]]$snpinfo))) == FALSE)
      hasrsq <- TRUE
  }

  snpinfo <- list()
  numsnpinfo <- 0
  snpinfonames <- character(0)
  if (hasaaf == TRUE) {
    numsnpinfo <- numsnpinfo + 1
    snpinfo[[numsnpinfo]] <- matrix(rep(as.numeric(NA), nrow(mergedinfo$snps) * numgroups),
                                    nrow(mergedinfo$snps), numgroups)
    snpinfonames <- c(snpinfonames, "aaf")
  }
  if (hasmaf == TRUE) {
    numsnpinfo <- numsnpinfo + 1
    snpinfo[[numsnpinfo]] <- matrix(rep(as.numeric(NA), nrow(mergedinfo$snps) * numgroups),
                                    nrow(mergedinfo$snps), numgroups)
    snpinfonames <- c(snpinfonames, "maf")
  }
  if (hasavgcall == TRUE) {
    numsnpinfo <- numsnpinfo + 1
    snpinfo[[numsnpinfo]] <- matrix(rep(as.numeric(NA), nrow(mergedinfo$snps) * numgroups),
                                    nrow(mergedinfo$snps), numgroups)
    snpinfonames <- c(snpinfonames, "avgcall")
  }
  if (hasrsq == TRUE) {
    numsnpinfo <- numsnpinfo + 1
    snpinfo[[numsnpinfo]] <- matrix(rep(as.numeric(NA), nrow(mergedinfo$snps) * numgroups),
                                    nrow(mergedinfo$snps), numgroups)
    snpinfonames <- c(snpinfonames, "rsq")
  }
  if(numsnpinfo > 0) {
    names(snpinfo) <- snpinfonames
    currentgroup <- 1L
    for (i in 1:length(geneticinfo)) {
      setgroups <- 1L
      if (class(geneticinfo[[i]]$additionalinfo) == "bdose-info")
        setgroups <- geneticinfo[[i]]$additionalinfo$numgroups
      for (j in 1:numsnpinfo) {
        if (is.na(match(names(snpinfo)[j], names(geneticinfo[[i]]$snpinfo))) == FALSE) {
          snpinfo[[j]][snpsbtom[[i]], currentgroup:(currentgroup + setgroups - 1)] <-
            geneticinfo[[i]]$snpinfo[[names(snpinfo)[j]]][snpsmtob[[i]],]
        }
      }
      currentgroup <- currentgroup + setgroups
    }
  }

  return(snpinfo)
}

mergegeneticinfo <- function(mergefile,
                             geneticinfo,
                             allsnps) {
  usesfid <- geneticinfo[[1]]$usesfid
  samples <- geneticinfo[[1]]$samples
  snpidformat <- geneticinfo[[1]]$snpidformat
  snps <- geneticinfo[[1]]$snps

  for (i in 2:length(geneticinfo)) {
    if (geneticinfo[[i]]$usesfid != usesfid)
      stop("Some files use FID and others do not")
    if (snpidformat != geneticinfo[[i]]$snpidformat)
      snpidformat <- 0L
    samples <- rbind(samples, geneticinfo[[i]]$samples)
    snps <- merge(snps, geneticinfo[[i]]$snps, all = allsnps)
  }

  if (nrow(unique(samples)) != nrow(samples))
    stop("There are duplicate samples in the files to merge")

  chr1 <- snps$chromosome[1]
  onechr <- all(snps$chromosome == chr1)

  return (list(filename = mergefile,
               usesfid = usesfid,
               samples = samples,
               onechr = FALSE,
               snpidformat = snpidformat,
               snps = snps,
               snpinfo = list(),
               datasize = integer(0),
               indices = numeric(0)))
}

mergebdinfo <- function(bdinfo, format, subformat, onegroup) {
  if (onegroup) {
    numgroups <- 1L
    groups <- 0L
    for (i in 1:length(bdinfo))
      groups <- groups + nrow(bdinfo[[i]]$samples)
  } else {
    groups <- integer()
    for (i in 1:length(bdinfo))
      groups <- c(groups, bdinfo[[i]]$additionalinfo$groups)
    numgroups <- length(groups)
  }
  if (subformat == 0L) {
    subformat <- 2L
    for (i in 1:length(bdinfo)) {
      if (bdinfo[[i]]$additionalinfo$subformat == 1)
        subformat <- 1L
    }
  }
  retval <- list(format = format,
                 subformat = subformat,
                 headersize = 0,
                 numgroups = numgroups,
                 groups = groups)
  class(retval) <- "bdose-info"

  return (retval)
}

bdmerge <- function(mergefiles,
                    format = 4,
                    subformat = 0L,
                    bdfiles,
                    famfiles = character(),
                    mapfiles = character(),
                    onegroup = TRUE,
                    bdoptions = character(),
                    snpjoin = "inner") {
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

  if (missing(mergefiles) == TRUE)
    stop("No output files specified")
  if (is.character(mergefiles) == FALSE)
    stop("Output file names must be a character values")
  if (format == 4 & length(mergefiles) != 1)
    stop("Only one file name is needed when using format 4")
  if (format < 4 & length(mergefiles) != 3)
    stop("Three file names are required when using formats 1, 2, and 3")
  if (is.na(match("", mergefiles)) == FALSE)
    stop("Output file names cannot be blank")

  if (missing(bdfiles) == TRUE)
    stop("No files specified")
  if (is.character(bdfiles) == FALSE)
    stop("bdfiles must be a character vector")
  if (length(bdfiles) < 2)
    stop("At least two binary dosage files must be specified")

  if (is.character(famfiles) == FALSE)
    stop("famfiles must be a character vector")
  if (is.character(mapfiles) == FALSE)
    stop("mapfiles must be a character vector")
  if (length(famfiles) != 0 & length(mapfiles) == 0)
    stop("If famfiles is specified, mapfiles must be specified")
  if (length(famfiles) != 0 & length(mapfiles) == 0)
    stop("If mapfiles is specified, famfiles must be specified")
  if (length(famfiles) != 0 & (length(famfiles) != length(bdfiles) | length(mapfiles) != length(bdfiles)))
    stop("If famfiles and mapfiles are specified they must have the same length as bdfiles")

  if (is.logical(onegroup) == FALSE)
    stop("onegroup must be logical value")
  if (length(onegroup) != 1)
    stop("onegroups must be a logical vector of length 1")

  if (is.character(bdoptions) == FALSE)
    stop("bdoptions must be a character vector")
  if (onegroup == FALSE & length(bdoptions) > 0)
    stop("bdoptions can only be used if onegroup is TRUE")

  if (snpjoin != "inner" & snpjoin != "outer")
    stop("snpjoin must have a value of either \"inner\" or \"outer\"")

  allsnps <- TRUE
  if (snpjoin == "inner")
    allsnps <- FALSE

  bdinfo <- vector("list", length(bdfiles))
  for (i in 1:length(bdfiles)) {
    if (length(famfiles) > 0)
      bdinfo[[i]] <- getbdinfo(c(bdfiles[i], famfiles[i], mapfiles[i]))
    else
      bdinfo[[i]] <- getbdinfo(bdfiles[i])
  }
  # Testing - needs to be removed
  ######################################################
  bdinfo[[1]]$snps$snpid[3] <- bdinfo[[2]]$snps$snpid[1]
  bdinfo[[1]]$snps$snpid[4] <- bdinfo[[2]]$snps$snpid[3]
  bdinfo[[1]]$snps$snpid[1] <- bdinfo[[2]]$snps$snpid[4]
  bdinfo[[1]]$snps[7,] <- bdinfo[[2]]$snps[9,]
  bdinfo[[1]]$snps[9,] <- bdinfo[[2]]$snps[7,]
  ######################################################

  bdmergedinfo <- mergebdinfo(bdinfo = bdinfo,
                              format = format,
                              subformat = subformat,
                              onegroup = onegroup)

  mergedgeneticinfo <- mergegeneticinfo(mergefile = mergefiles[1],
                                        geneticinfo = bdinfo,
                                        allsnps = allsnps)
  mergedgeneticinfo$additionalinfo <- bdmergedinfo

  snpsbtom <- list()
  snpsmtob <- list()
  snpsbtomNA <- list()
  snpsmtobNA <- list()
  for (i in 1:length(bdinfo)) {
    snpsbtomNA[[i]] <- prodlim::row.match(bdinfo[[i]]$snps, mergedgeneticinfo$snps)
    snpsbtom[[i]] <- snpsbtomNA[[i]][!is.na(snpsbtomNA[[i]])]
    snpsmtobNA[[i]] <- prodlim::row.match(mergedgeneticinfo$snps, bdinfo[[i]]$snps)
    snpsmtob[[i]] <- snpsmtobNA[[i]][!is.na(snpsmtobNA[[i]])]
  }

  mergedgeneticinfo$snpinfo <- mergesnpinfo(mergedinfo = mergedgeneticinfo,
                                            numgroups = bdmergedinfo$numgroups,
                                            geneticinfo = bdinfo,
                                            snpsbtom = snpsbtom,
                                            snpsmtob = snpsmtob)

  WriteBinaryDosageHeader(format = bdmergedinfo$format,
                          subformat = bdmergedinfo$subformat,
                          filename = mergefiles,
                          genefileinfo = mergedgeneticinfo,
                          bdoptions = bdoptions)
  headerinfo <- ReadBinaryDosageHeader(filename = mergefiles)
  bdwriteinfo <- AllocateBinaryDosageWriteMemory(headerinfo = headerinfo)

  dosage <- numeric(nrow(mergedgeneticinfo$samples))
  p0 <- numeric(nrow(mergedgeneticinfo$samples))
  p1 <- numeric(nrow(mergedgeneticinfo$samples))
  p2 <- numeric(nrow(mergedgeneticinfo$samples))
  us <- integer(2 * nrow(mergedgeneticinfo$samples))
  dosaget <- numeric(nrow(mergedgeneticinfo$samples))
  p0t <- numeric(nrow(mergedgeneticinfo$samples))
  p1t <- numeric(nrow(mergedgeneticinfo$samples))
  p2t <- numeric(nrow(mergedgeneticinfo$samples))
  ust <- integer(2 * nrow(mergedgeneticinfo$samples))
  startgroup <- integer(length(bdinfo))
  endgroup <- integer(length(bdinfo))
  startgroup[1] <- 1
  endgroup[1] <- nrow(bdinfo[[1]]$samples)
  for (i in 2:length(bdinfo)) {
    startgroup[i] <- endgroup[i - 1] + 1
    endgroup[i] <- endgroup[i - 1] + nrow(bdinfo[[i]]$samples)
  }

  for (i in 1:nrow(mergedgeneticinfo$snps)) {
#  for (i in 1:1) {
    for (j in 1:length(bdinfo)) {
      dosage[1:nrow(mergedgeneticinfo$samples)] <- NA
      p0[1:nrow(mergedgeneticinfo$samples)] <- NA
      p1[1:nrow(mergedgeneticinfo$samples)] <- NA
      p2[1:nrow(mergedgeneticinfo$samples)] <- NA
      if (is.na(snpsmtobNA[[j]][i]) == FALSE)
        ReadBinaryDosageData(bdinfo[[j]], snpsmtobNA[[j]][i], dosage, p0, p1, p2, us)
      dosaget[startgroup[j]:endgroup[j]] <- dosage[1:nrow(bdinfo[[j]]$samples)]
      p0t[startgroup[j]:endgroup[j]] <- p0[1:nrow(bdinfo[[j]]$samples)]
      p1t[startgroup[j]:endgroup[j]] <- p1[1:nrow(bdinfo[[j]]$samples)]
      p2t[startgroup[j]:endgroup[j]] <- p2[1:nrow(bdinfo[[j]]$samples)]
    }
    WriteBinaryDosageData(dosaget, p0t, p1t, p2t, bdwriteinfo)
    print(mean(dosaget, na.rm = TRUE) / 2.)
  }
  WriteBinaryDosageIndices(writeinfo = bdwriteinfo)

  return (list(bdmergedinfo = bdmergedinfo,
               mergedgeneticinfo = mergedgeneticinfo,
               snpsmtob = snpsmtob,
               snpsbtom = snpsbtom,
               snpsmtobNA = snpsmtobNA,
               snpsbtomNA = snpsbtomNA,
               dosage = dosage,
               p0 = p0,
               p1 = p1,
               p2 = p2,
               dosaget = dosaget,
               p0t = p0t,
               p1t = p1t,
               p2t = p2t,
               bdinfo = bdinfo))
}
