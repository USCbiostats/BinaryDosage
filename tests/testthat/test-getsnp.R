test_that("getsnp", {
  vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcfinfo <- getvcfinfo(vcffiles = vcf1afile)
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdinfo <- getbdinfo(bdfiles = vcf1abdfile)

  expect_error(getsnp(snp = 1),
               "bdinfo missing")
  expect_error(getsnp(bdinfo = 1,
                      snp = 1),
               "bdinfo is not of class genetic-info")
  expect_error(getsnp(bdinfo = vcfinfo,
                      snp = 1),
               "bdinfo does not contain information about a binary dosage file")

  expect_error(getsnp(bdinfo = bdinfo),
               "No SNP specified")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1:3),
               "snp must be of length one")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1.3),
               "snp must be a character or integer value")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = TRUE),
               "snp must be a character or integer value")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = "xyz"),
               "Cannot find SNP in bdinfo")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 100),
               "snp value out or range")

  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1,
                      FALSE),
               NA)
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = "1:11000:T:C",
                      TRUE),
               NA)
})

test_that("getsnp format5", {
  skip_if_not_installed("vcfppR")

  vcfgzfile <- system.file("extdata", "set1a.vcf.gz", package = "BinaryDosage")
  bdose_file  <- tempfile(fileext = ".bdose")
  bdinfo_file <- tempfile(fileext = ".bdinfo")

  expect_error(vcftobd(vcffile = vcfgzfile,
                       bdose_file  = bdose_file,
                       bdinfo_file = bdinfo_file),
               NA)

  bd5info <- getbd5info(bdose_file = bdose_file, bdinfo_file = bdinfo_file)
  n_snps <- nrow(bd5info$snps)
  expect_equal(n_snps, 10L)

  for (i in seq_len(n_snps)) {
    ref    <- getbd5snp(bd5info, i)
    result <- getsnp(bdinfo = bd5info, snp = i, dosageonly = FALSE)

    expect_equal(result$dosage, ref$dosage, tolerance = 5e-5,
                 label = paste0("dosage SNP ", i))
    expect_equal(result$p0, ref$p0, tolerance = 5e-5,
                 label = paste0("p0 SNP ", i))
    expect_equal(result$p1, ref$p1, tolerance = 5e-5,
                 label = paste0("p1 SNP ", i))
    expect_equal(result$p2, ref$p2, tolerance = 5e-5,
                 label = paste0("p2 SNP ", i))
  }
})
