test_that("merge", {
  expect_error(bdmerge(bdfiles = c("file1", "file2")),
               "No output files specified")
  expect_error(bdmerge(mergefiles = 1L,
                       bdfiles = c("file1", "file2")),
               "Output file names must be a character values")
  expect_error(bdmerge(mergefiles = c("file1", "file2", "file3"),
                       bdfiles = c("file1", "file2")),
               "Only one file name is needed when using format 4")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 3,
                       bdfiles = c("file1", "file2")),
               "Three file names are required when using formats 1, 2, and 3")
  expect_error(bdmerge(mergefiles = "",
                       bdfiles = c("file1", "file2")),
               "Output file names cannot be blank")

  expect_error(bdmerge(mergefiles = "file1"),
               "No files specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = 1),
               "bdfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = "file1"),
               "At least two binary dosage files must be specified")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = 1,
                       mapfiles = c("file1", "file2")),
               "famfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2"),
                       mapfiles = 1),
               "mapfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2")),
               "If famfiles is specified, mapfiles must be specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       mapfiles = c("file1", "file2")),
               "If mapfiles is specified, famfiles must be specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2"),
                       mapfiles = "file1"),
               "If famfiles and mapfiles are specified they must have the same length as bdfiles")

  expect_error(bdmerge(mergefiles = "file1",
                       format = "a",
                       bdfiles = c("file1", "file2")),
               "format must be an integer value")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 1:2,
                       bdfiles = c("file1", "file2")),
               "format must be an integer vector of length 1")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 1.2,
                       bdfiles = c("file1", "file2")),
               "format must be an integer")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 5,
                       bdfiles = c("file1", "file2")),
               "format must be an integer value from 1 to 4")

  expect_error(bdmerge(mergefiles = "file1",
                       subformat = "a",
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 1:2,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer vector of length 1")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 1.2,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 5,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value from 0 to 4")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 2,
                       subformat = 3,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value from 0 to 2 for formats 1 and 2")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = 1),
               "onegroup must be logical value")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = c(TRUE, TRUE)),
               "onegroup must be a logical vector of length 1")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       bdoptions = 1),
               "bdoptions must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = FALSE,
                       bdoptions = "aaf"),
               "bdoptions can only be used if onegroup is TRUE")
  expect_error(bdmerge(mergefiles = c("file1", "file2", "file3"),
                       bdfiles = c("file1", "file2"),
                       format = 3,
                       bdoptions = "aaf"),
               "bdoptions can only be used if format = 4")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       bdoptions = "abc"),
               "Only valid bdoptions are aaf, maf, and rsq")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       snpjoin = "abc"),
               "snpjoin must have a value of either \"inner\" or \"outer\"")

  bdvcf1afile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdvcf1bfile <- system.file("extdata", "vcf1b.bdose", package = "BinaryDosage")
  mergefiles <- tempfile()

  expect_error(BinaryDosage:::bdmerge(mergefiles = mergefiles,
                                      bdfiles = c(bdvcf1afile, bdvcf1bfile),
                                      bdoptions = "maf"),
               NA)
  expect_error(bdinfo <- getbdinfo(mergefiles), NA)
})
