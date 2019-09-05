mergegeneticinfo <- function(geneticinfo, allsnps) {
  usesfid <- geneticinfo[[1]]$usesfid
  for (i in 2:length(geneticinfo)) {
    if (geneticinfo[[i]]$usesfid != usesfid)
      stop("Some files use FID and others do not")
  }

  samples <- geneticinfo[[1]]$samples
  for (i in 2:length(bdinfo))
    samples <- rbind(samples, geneticinfo[[i]]$samples)
#  print(nrow(unique(samples)))
#  print(nrow(samples))
  if (nrow(unique(samples)) != nrow(samples))
    print("There are duplicate samples in the merged file")

  snpidformat <- geneticinfo[[1]]$snpidformat
  for (i in 2:length(geneticinfo)) {
    if (geneticinfo[[i]]$snpidformat != snpidformat)
      stop("All files must use same SNP ID format")
  }

  snps <- geneticinfo[[1]]$snps
  snps$x <- 1:nrow(snps)
  colnames(snps)[6] <- "file1"
  for (i in 2:length(geneticinfo)) {
    mergesnps <- geneticinfo[[1]]$snps
    mergesnps$x <- 1:nrow(mergesnps)
    colnames(mergesnps)[6] <- paste("file", i, sep = "")
    snps <- merge(x = snps,
                  y = mergesnps,
                  all = allsnps)
  }
  snps <- snps[order(snps$location),]
  if (anyNA(snps[,6:ncol(snps)]) == TRUE)
    stop("SNPs are not the same in all files")
  onechr <- FALSE
  if (length(unique(snps$chromosome)) == 1)
    onechr <- TRUE

  snpinfo <- geneticinfo[[1]]$snpinfo
  print(length(snpinfo))
  for (i in 2:length(geneticinfo)) {
    print(length(geneticinfo[[i]]$snpinfo))
    if (length(snpinfo) != length(geneticinfo[[i]]$snpinfo))
      stop("All files must have the same SNP information")
    if (length(snpinfo) > 0){
      if (colnames(geneticinfo[[i]]$snpinfo) != colnames(snpinfo))
        stop("All files must have the same SNP information")
      for (j in 1:length(snpinfo))
        snpinfo[[j]] <- c(snpinfo[[j]], geneticinfo[[i]]$snpinfo[[j]])
    }
  }

  return (list(filename = "",
               usesfid = usesfid,
               samples = samples,
               onechr = onechr,
               snpidformat = snpidformat,
               snps = snps,
               snpinfo = snpinfo,
               datasize = integer(0),
               indices = numeric(0)))
}

mergebdinfo <- function(bdinfo) {

}

bdmerge <- function(bdfiles,
                    famfiles = character(),
                    mapfiles = character(),
                    onegroup = TRUE,
                    snpjoin = "inner") {

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

  return (bdinfo)
}
