bdmerge <- function(bdfiles) {
  print(length(bdfiles))
  for (i in 1:length(bdfiles))
    print(bdfiles[[i]])
  bdinfo <- vector("list", length(bdfiles))
  for (i in 1:length(bdfiles))
    bdinfo[[i]] <- getbdinfo(bdfiles[[i]])
  return (bdinfo)
}
