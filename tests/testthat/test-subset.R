## Helper: compute MAF of a SNP over a subset of subject indices
snp_maf <- function(bdinfo, snp_idx, samp_idx = NULL) {
  snp <- getsnp(bdinfo, snp_idx, dosageonly = TRUE)
  ds  <- if (is.null(samp_idx)) snp$dosage else snp$dosage[samp_idx]
  aaf <- mean(ds, na.rm = TRUE) / 2.0
  min(aaf, 1.0 - aaf)
}

test_that("subsetbd errors", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")

  expect_error(subsetbd(bdose_file = "x.bdose", bdinfo_file = "x.bdinfo"),
               "No input binary dosage files specified")
  expect_error(subsetbd(bdfiles = bdfile, bdinfo_file = "x.bdinfo"),
               "No output .bdose file specified")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose"),
               "No output .bdinfo file specified")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo"),
               "At least one of")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        locations = c(10000, 11000),
                        startloc = 10000),
               "locations and startloc/endloc cannot both be specified")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        startloc = 10000),
               "endloc must be specified")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        endloc = 19000),
               "startloc must be specified")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        startloc = 15000, endloc = 10000),
               "startloc must be less than or equal to endloc")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        minmaf = "high"),
               "minmaf must be a single numeric value")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        minmaf = 0.6),
               "minmaf must be between 0 and 0.5")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        subjectids = c("X1", "X2")),
               "No subjects match")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = "x.bdose",
                        bdinfo_file = "x.bdinfo",
                        minmaf = 0.5),
               "No SNPs pass")
})

test_that("subsetbd subject filter", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  old_info <- getbdinfo(bdfile)
  keep_ids <- old_info$samples$sid[1:30]

  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = bdose,
                        bdinfo_file = bdinf, subjectids = keep_ids), NA)

  new_info <- getbd5info(bdose, bdinf)
  expect_equal(nrow(new_info$samples), 30L)
  expect_equal(new_info$samples$sid, keep_ids)
  expect_equal(nrow(new_info$snps),   nrow(old_info$snps))

  # Verify dosage values match for a sample of SNPs
  for (i in c(1L, 5L, 10L)) {
    old_snp <- getsnp(old_info, i, dosageonly = TRUE)
    new_snp <- getsnp(new_info, i, dosageonly = TRUE)
    expect_equal(new_snp$dosage, old_snp$dosage[1:30], tolerance = 5e-5,
                 label = paste0("subject filter dosage SNP ", i))
  }
})

test_that("subsetbd location vector filter", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  old_info <- getbdinfo(bdfile)
  keep_locs <- c(12000, 15000, 17000)   # SNPs 3, 6, 8

  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = bdose,
                        bdinfo_file = bdinf, locations = keep_locs), NA)

  new_info <- getbd5info(bdose, bdinf)
  expect_equal(nrow(new_info$snps), length(keep_locs))
  expect_equal(new_info$snps$location, keep_locs)
  expect_equal(nrow(new_info$samples), nrow(old_info$samples))

  # Verify data for each kept SNP
  orig_idx <- match(keep_locs, old_info$snps$location)
  for (j in seq_along(keep_locs)) {
    old_snp <- getsnp(old_info, orig_idx[j], dosageonly = TRUE)
    new_snp <- getsnp(new_info, j, dosageonly = TRUE)
    expect_equal(new_snp$dosage, old_snp$dosage, tolerance = 5e-5,
                 label = paste0("location filter dosage SNP loc=", keep_locs[j]))
  }
})

test_that("subsetbd location range filter", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  old_info <- getbdinfo(bdfile)
  # startloc=13000, endloc=16000 keeps SNPs 4,5,6,7 (locs 13000-16000)
  exp_locs <- c(13000, 14000, 15000, 16000)

  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = bdose,
                        bdinfo_file = bdinf,
                        startloc = 13000, endloc = 16000), NA)

  new_info <- getbd5info(bdose, bdinf)
  expect_equal(new_info$snps$location, exp_locs)

  orig_idx <- match(exp_locs, old_info$snps$location)
  for (j in seq_along(exp_locs)) {
    old_snp <- getsnp(old_info, orig_idx[j], dosageonly = TRUE)
    new_snp <- getsnp(new_info, j,            dosageonly = TRUE)
    expect_equal(new_snp$dosage, old_snp$dosage, tolerance = 5e-5,
                 label = paste0("range filter dosage loc=", exp_locs[j]))
  }
})

test_that("subsetbd MAF filter", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  old_info <- getbdinfo(bdfile)
  min_maf  <- 0.25

  # Determine expected kept SNPs
  keep_snps <- which(vapply(seq_len(nrow(old_info$snps)),
                            function(i) snp_maf(old_info, i) >= min_maf,
                            logical(1L)))

  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles = bdfile, bdose_file = bdose,
                        bdinfo_file = bdinf, minmaf = min_maf), NA)

  new_info <- getbd5info(bdose, bdinf)
  expect_equal(nrow(new_info$snps), length(keep_snps))
  expect_equal(new_info$snps$location,
               old_info$snps$location[keep_snps])

  for (j in seq_along(keep_snps)) {
    old_snp <- getsnp(old_info, keep_snps[j], dosageonly = TRUE)
    new_snp <- getsnp(new_info, j,             dosageonly = TRUE)
    expect_equal(new_snp$dosage, old_snp$dosage, tolerance = 5e-5,
                 label = paste0("maf filter dosage SNP ", keep_snps[j]))
  }
})

test_that("subsetbd combined filters", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  old_info <- getbdinfo(bdfile)

  keep_ids <- old_info$samples$sid[1:30]
  samp_idx <- 1:30
  min_maf  <- 0.25

  # Expected: range 13000-18000, maf >= 0.25 over first 30 subjects
  exp_locs <- c(13000, 14000, 15000, 16000, 17000, 18000)
  loc_snps <- match(exp_locs, old_info$snps$location)
  keep_snps <- loc_snps[vapply(loc_snps,
                               function(i) snp_maf(old_info, i, samp_idx) >= min_maf,
                               logical(1L))]

  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles     = bdfile,
                        bdose_file  = bdose,
                        bdinfo_file = bdinf,
                        minmaf      = min_maf,
                        startloc    = 13000,
                        endloc      = 18000,
                        subjectids  = keep_ids), NA)

  new_info <- getbd5info(bdose, bdinf)
  expect_equal(nrow(new_info$samples), 30L)
  expect_equal(new_info$samples$sid, keep_ids)
  expect_equal(new_info$snps$location, old_info$snps$location[keep_snps])

  for (j in seq_along(keep_snps)) {
    old_snp <- getsnp(old_info, keep_snps[j], dosageonly = TRUE)
    new_snp <- getsnp(new_info, j,             dosageonly = TRUE)
    expect_equal(new_snp$dosage, old_snp$dosage[samp_idx], tolerance = 5e-5,
                 label = paste0("combined filter dosage SNP ", keep_snps[j]))
  }
})

test_that("subsetbd format 5 input", {
  # Verify subsetbd also accepts format 5 as input
  bdfile4   <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdose5_in <- tempfile(fileext = ".bdose")
  bdinf5_in <- tempfile(fileext = ".bdinfo")
  updatebd(bdfiles = bdfile4, bdose_file = bdose5_in, bdinfo_file = bdinf5_in)

  bdose_out <- tempfile(fileext = ".bdose")
  bdinf_out <- tempfile(fileext = ".bdinfo")
  expect_error(subsetbd(bdfiles     = c(bdose5_in, bdinf5_in),
                        bdose_file  = bdose_out,
                        bdinfo_file = bdinf_out,
                        startloc    = 12000,
                        endloc      = 14000), NA)

  new_info <- getbd5info(bdose_out, bdinf_out)
  expect_equal(new_info$snps$location, c(12000, 13000, 14000))
})
