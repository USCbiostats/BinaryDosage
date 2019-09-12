test_that("getbdinfo", {
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")

  expect_error(getbdinfo(),
               "No binary dosage files specified")
  expect_error(getbdinfo(bdfiles = 1),
               "bdfiles must be a character vector")
  expect_error(getbdinfo(bdfiles = character()),
               "bdfiles must be a character vector of length 1 or 3")
  expect_error(getbdinfo(bdfiles = ""),
               "No binary dosage file specified")
  expect_error(getbdinfo(bdfiles = c("file1", "file2", "")),
               "bdfiles contains empty strings")
  expect_error(getbdinfo(bdfiles = c(vcf1abdfile, "file2", "file3")),
               "Binary dosage file format 4 does not use family and map files")
})

test_that("getvcfinfo", {
  expect_error(getvcfinfo(),
               "No VCF file specified")
  expect_error(getvcfinfo(vcffiles = 1),
               "vcfiles must be a character value")
  expect_error(getvcfinfo(vcffiles = c("file1", "file2", "file3")),
               "vcffiles must be a character vector of length 1 or 2")
  expect_error(getvcfinfo(vcffiles = character()),
               "vcffiles must be a character vector of length 1 or 2")
  expect_error(getvcfinfo(vcffiles = c("", "file2")),
               "No VCF file specified")

  expect_error(getvcfinfo(vcffiles = "file1",
                          gz = 1L),
               "gz must be a logical value")
  expect_error(getvcfinfo(vcffiles = "file1",
                          gz = c(TRUE, TRUE)),
               "gz must be a logical vector of length 1")
  expect_error(getvcfinfo(vcffiles = "file1",
                          index = 1L),
               "index must be a logical value")
  expect_error(getvcfinfo(vcffiles = "file1",
                          index = c(TRUE, TRUE)),
               "index must be a logical vector of length 1")
  expect_error(getvcfinfo(vcffiles = "file1",
                          gz = TRUE,
                          index = TRUE),
               "Indexing gzipped files is not supported.")

  expect_error(getvcfinfo(vcffiles = "file1",
                          snpidformat = ""),
               "snpidformat must be an intger value")
  expect_error(getvcfinfo(vcffiles = "file1",
                          snpidformat = c(1L, 2L)),
               "snpidformat must be a singe integer value")
  expect_error(getvcfinfo(vcffiles = "file1",
                          snpidformat = 1.1),
               "snpidformat must be an integer value")
  expect_error(getvcfinfo(vcffiles = "file1",
                          snpidformat = 3),
               "snpidformat must have a value of 0, 1, or 2")
})

test_that("getgeninfo", {
  expect_error(getgeninfo(),
               "No gen file specified")
  expect_error(getgeninfo(genfiles = 1L),
               "genfiles must be a character value")
  expect_error(getgeninfo(genfiles = c("file1", "file2", "file3")),
               "genfiles must be a character vector of length 1 or 2")
  expect_error(getgeninfo(genfiles = character()),
               "genfiles must be a character vector of length 1 or 2")
  expect_error(getgeninfo(genfiles = c("", "file2")),
               "No gen file specified")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpcolumns = ""),
               "snpcolumns must be an integer vector")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpcolumns = 1L:3L),
               "snpcolumns must be an integer vector of length 5")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpcolumns = c(1L, 0L, 2L, 3L, 4L)),
               "snpcolumns values other than chromosome must be positive integers")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpcolumns = c(-2L, 1L, 2L, 3L, 4L)),
               "snpcolumns chromosome value must be -1, or a non-negative integer")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          startcolumn = ""),
               "startcolumn must be an integer value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          startcolumn = 1L:3L),
               "startcolumn must be an integer vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          startcolumn = -1L),
               "startcolumn must be a positive integer")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          startcolumn = 1L),
               "startcolumn value must be larger than any value in snpcolumns")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          impformat = ""),
               "impformat must be an integer value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          impformat = 1L:3L),
               "impformat must be an integer vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          impformat = 4L),
               "impformat must have a value of 1, 2, or 3")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          chromosome = 1),
               "chromosome must be a character variable")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          chromosome = c("a", "b")),
               "chromosome must be a character vector of length 0 or 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpcolumns = c(-1L,2L:5L)),
               "No chromosome column or value provided")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          chromosome = "1"),
               "Both chromosome column and chromosome value provided")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          header = 1),
               "header must be a logical value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          header = c(TRUE, TRUE)),
               "header must be a logical vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", ""),
                          header = FALSE),
               "File has no header and no sample file is provided")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          header = TRUE),
               "header = TRUE and a sample file is provided")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          gz = 1L),
               "gz must be a logical value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          gz = c(TRUE, TRUE)),
               "gz must be a logical vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          index = 1L),
               "index must be a logical value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          index = c(TRUE, TRUE)),
               "index must be a logical vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          gz = TRUE,
                          index = TRUE),
               "Indexing gzipped files is not supported.")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpidformat = ""),
               "snpidformat must be an integer value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpidformat = 1L:3L),
               "snpidformat must be an integer vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          snpidformat = 5L),
               "snpidformat must have a value of 0, 1, 2, or 3")

  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          sep = 1L),
               "sep must be a character value")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          sep = c("a", "b")),
               "sep must be a character vector of length 1")
  expect_error(getgeninfo(genfiles = c("file1", "file2"),
                          sep = ""),
               "no value of sep provided")

})
