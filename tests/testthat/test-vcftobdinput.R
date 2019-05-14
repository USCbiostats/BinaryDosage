context("test-vcftobdinput")

test_that("VCFtoBD input valid", {
  vcfFile <- "junk.vcf"
  bdFile <- "junk.bdose"

  expect_error(VCFtoBD(bdFile = bdFile),
               "No VCF file specified")
  expect_error(VCFtoBD(vcfFile = 1,
                       bdFile = bdFile),
               "vcfFile must be a character value")
  expect_error(VCFtoBD(vcfFile = c("a", "b"),
                       bdFile = bdFile),
               "vcfFile must be a single character value")

  expect_error(VCFtoBD(vcfFile = vcfFile),
               "No output file specified")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = 1),
               "bdFile must be a character value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = c("a", "b")),
               "bdFile must be a single character value")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       format = ""),
               "format must be an integer value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       format = c(1,2)),
               "format must be a single integer value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       format = 1.2),
               "format must be an integer")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       format = 5),
               "format must be an integer value from 1 to 4")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       version = ""),
               "version must be an integer value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       version = c(1,2)),
               "version must be a single integer value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       version = 1.2),
               "version must be an integer")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       version = 5),
               "version must be an integer value from 0 to 2")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = 1,
                       mapFile = "map",
                       format = 3),
               "famFile must be a character value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = c("a", "b"),
                       mapFile = "map",
                       format = 3),
               "famFile must be a single character value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "",
                       mapFile = "map",
                       format = 3),
               "famFile must be specified for format 3")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "fam",
                       mapFile = 1,
                       format = 2),
               "mapFile must be a character value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "fam",
                       mapFile = c("a", "b"),
                       format = 2),
               "mapFile must be a single character value")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "fam",
                       mapFile = "",
                       format = 2),
               "mapFile must be specified for format 2")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "fam",
                       mapFile = "",
                       format = 4),
               "Value or famFile specified for format 4")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = 1,
                       mapFile = "",
                       format = 4),
               "Value or famFile specified for format 4")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = c("a", "b"),
                       mapFile = "",
                       format = 4),
               "Value or famFile specified for format 4")

  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "",
                       mapFile = "map",
                       format = 4),
               "Value or mapFile specified for format 4")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "",
                       mapFile = 1,
                       format = 4),
               "Value or mapFile specified for format 4")
  expect_error(VCFtoBD(vcfFile = vcfFile,
                       bdFile = bdFile,
                       famFile = "",
                       mapFile = c("a", "b"),
                       format = 4),
               "Value or mapFile specified for format 4")
})
