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
