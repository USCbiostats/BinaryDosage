test_that("vcftobd", {
  expect_error(vcftobd(bdfiles = "file1"),
               "No VCF file specified")
  expect_error(vcftobd(vcffiles = "file1"),
               "No output files specified")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = 1),
               "bdfiles must be a vector of characters")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = ""),
               "format must be an integer value")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 1:2),
               "format must be an integer vector of length 1")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 1.2),
               "format must be an integer")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 5),
               "format must be an integer value from 1 to 4")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = ""),
               "subformat must be an integer value")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 1:2),
               "subformat must be an integer vector of length 1")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 1.2),
               "subformat must be an integer")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 5),
               "subformat must be an integer value from 0 to 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 2,
                       subformat = 3),
               "subformat must be an integer value from 0 to 2 for formats 1 and 2")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3"),
                       format = 4),
               "Only one output file name is needed when using format 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3"),
                       format = 3),
               "Three output file names are required when using formats 1, 2, and 3")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3", ""),
                       format = 3),
               "Output file names cannot be blank")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       bdoptions = 4),
               "bdoptions must be a character array")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3", "file4"),
                       format = 3,
                       bdoptions = "aaf"),
               "bdoptions can only be used with format 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       bdoptions = "abc"),
               "Only valid bdoptions are aaf, maf, and rsq")

  vcf1afile = system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
  bdfiles <- tempfile()
  expect_error(vcftobd(vcffiles = c(vcf1afile, vcf1ainfo),
                       bdfiles = bdfiles,
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles = bdfiles), NA)

  vcf2afile = system.file("extdata", "set2a.vcf", package = "BinaryDosage")
  bdfiles2 <- tempfile()
  expect_error(vcftobd(vcffiles = c(vcf2afile),
                       bdfiles = bdfiles2),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles = bdfiles2), NA)

})

test_that("gentobd", {
  expect_error(gentobd(bdfiles = "file1"),
               "No gen file specified")
  expect_error(gentobd(genfiles = "file1"),
               "No output files specified")

  gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
  gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
  bdfiles3a <- tempfile()
  expect_error(gentobd(genfiles = c(gen3afile, gen3asample),
                       snpcolumns = c(0L, 2L:5L),
                       bdfiles = bdfiles3a,
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles3a), NA)

  gen2afile <- system.file("extdata", "set2a.imp", package = "BinaryDosage")
  gen2asample <- system.file("extdata", "set2a.sample", package = "BinaryDosage")
  bdfiles2a <- tempfile()
  expect_error(gentobd(genfiles = c(gen2afile, gen2asample),
                       snpcolumns = c(1L, 3L, 2L, 4L, 5L),
                       impformat = 1L,
                       bdfiles = bdfiles2a),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles2a), NA)})
