context("test-getvcfinfoinput")

test_that("Valid GetVCFInfoInput", {
  expect_error(GetVCFInfo(),
               "No VCF file specified")
  expect_error(GetVCFInfo(vcfFile = ""),
               "No VCF file specified")
  expect_error(GetVCFInfo(vcfFile = 1),
               "vcfFile must be a character value")
  expect_error(GetVCFInfo(vcfFile = c("a", "b")),
               "vcfFile must be a single character value")

  expect_error(GetVCFInfo(vcfFile = "junk.txt",
                          index = 1),
               "index must be a logical value")
  expect_error(GetVCFInfo(vcfFile = "junk.txt",
                          index = c(TRUE, TRUE)),
               "index must be a single logical value")
})
