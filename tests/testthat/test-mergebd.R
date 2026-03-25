## Helper: convert format-4 extdata file to Format 5
make_bd5 <- function(bdfile4) {
  bdose <- tempfile(fileext = ".bdose")
  bdinf <- tempfile(fileext = ".bdinfo")
  updatebd(bdfiles = bdfile4, bdose_file = bdose, bdinfo_file = bdinf)
  list(bdose = bdose, bdinf = bdinf, info = getbd5info(bdose, bdinf))
}

test_that("mergebd errors", {
  bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  f5     <- make_bd5(bdfile)

  expect_error(mergebd(bdinfo_files = c(f5$bdinf, f5$bdinf),
                       bdose_file = "x.bdose", bdinfo_file = "x.bdinfo"),
               "No input .bdose files specified")
  expect_error(mergebd(bdose_files  = c(f5$bdose, f5$bdose),
                       bdose_file = "x.bdose", bdinfo_file = "x.bdinfo"),
               "No input .bdinfo files specified")
  expect_error(mergebd(bdose_files  = c(f5$bdose, f5$bdose),
                       bdinfo_files = c(f5$bdinf, f5$bdinf),
                       bdinfo_file  = "x.bdinfo"),
               "No output .bdose file specified")
  expect_error(mergebd(bdose_files  = c(f5$bdose, f5$bdose),
                       bdinfo_files = c(f5$bdinf, f5$bdinf),
                       bdose_file   = "x.bdose"),
               "No output .bdinfo file specified")
  expect_error(mergebd(bdose_files  = c(f5$bdose, f5$bdose),
                       bdinfo_files = c(f5$bdinf, f5$bdinf, f5$bdinf),
                       bdose_file   = "x.bdose", bdinfo_file = "x.bdinfo"),
               "must have the same length")
  expect_error(mergebd(bdose_files  = f5$bdose,
                       bdinfo_files = f5$bdinf,
                       bdose_file   = "x.bdose", bdinfo_file = "x.bdinfo"),
               "At least two input files")

  # Both subjects and SNPs intersect (same file twice) -> error
  expect_error(mergebd(bdose_files  = c(f5$bdose, f5$bdose),
                       bdinfo_files = c(f5$bdinf, f5$bdinf),
                       bdose_file   = "x.bdose", bdinfo_file = "x.bdinfo"),
               "Cannot merge: subject IDs and SNP IDs both overlap")
})

test_that("mergebd subject merge (two files)", {
  bdfile   <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  f5       <- make_bd5(bdfile)
  all_info <- f5$info
  all_sids <- all_info$samples$sid   # I1..I60

  # Split into two non-overlapping subject sets (same SNPs)
  bdose_a <- tempfile(fileext = ".bdose"); bdinf_a <- tempfile(fileext = ".bdinfo")
  bdose_b <- tempfile(fileext = ".bdose"); bdinf_b <- tempfile(fileext = ".bdinfo")
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_a,
           bdinfo_file = bdinf_a, subjectids = all_sids[1:30])
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_b,
           bdinfo_file = bdinf_b, subjectids = all_sids[31:60])

  bdose_out <- tempfile(fileext = ".bdose"); bdinf_out <- tempfile(fileext = ".bdinfo")
  expect_error(mergebd(bdose_files  = c(bdose_a, bdose_b),
                       bdinfo_files = c(bdinf_a, bdinf_b),
                       bdose_file   = bdose_out,
                       bdinfo_file  = bdinf_out), NA)

  merged <- getbd5info(bdose_out, bdinf_out)

  # All 60 subjects, all 10 SNPs
  expect_equal(nrow(merged$samples), 60L)
  expect_equal(merged$samples$sid, all_sids)
  expect_equal(nrow(merged$snps),   nrow(all_info$snps))
  expect_equal(merged$snps$location, all_info$snps$location)

  # Verify dosage data for all SNPs
  for (i in seq_len(nrow(merged$snps))) {
    orig_snp   <- getsnp(all_info, i, dosageonly = TRUE)
    merged_snp <- getsnp(merged,   i, dosageonly = TRUE)
    expect_equal(merged_snp$dosage, orig_snp$dosage, tolerance = 5e-5,
                 label = paste0("subject merge dosage SNP ", i))
  }
})

test_that("mergebd subject merge (three files)", {
  bdfile   <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  f5       <- make_bd5(bdfile)
  all_sids <- f5$info$samples$sid

  bdose_a <- tempfile(fileext = ".bdose"); bdinf_a <- tempfile(fileext = ".bdinfo")
  bdose_b <- tempfile(fileext = ".bdose"); bdinf_b <- tempfile(fileext = ".bdinfo")
  bdose_c <- tempfile(fileext = ".bdose"); bdinf_c <- tempfile(fileext = ".bdinfo")
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_a,
           bdinfo_file = bdinf_a, subjectids = all_sids[1:20])
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_b,
           bdinfo_file = bdinf_b, subjectids = all_sids[21:40])
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_c,
           bdinfo_file = bdinf_c, subjectids = all_sids[41:60])

  bdose_out <- tempfile(fileext = ".bdose"); bdinf_out <- tempfile(fileext = ".bdinfo")
  expect_error(mergebd(bdose_files  = c(bdose_a, bdose_b, bdose_c),
                       bdinfo_files = c(bdinf_a, bdinf_b, bdinf_c),
                       bdose_file   = bdose_out,
                       bdinfo_file  = bdinf_out), NA)

  merged <- getbd5info(bdose_out, bdinf_out)
  expect_equal(nrow(merged$samples), 60L)
  expect_equal(merged$samples$sid, all_sids)
  expect_equal(nrow(merged$snps),  nrow(f5$info$snps))
})

test_that("mergebd SNP merge (two files)", {
  bdfile   <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  f5       <- make_bd5(bdfile)
  all_info <- f5$info

  # Split into two non-overlapping SNP sets (same subjects)
  locs     <- all_info$snps$location   # 10000..19000
  bdose_a <- tempfile(fileext = ".bdose"); bdinf_a <- tempfile(fileext = ".bdinfo")
  bdose_b <- tempfile(fileext = ".bdose"); bdinf_b <- tempfile(fileext = ".bdinfo")
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_a,
           bdinfo_file = bdinf_a, locations = locs[1:5])
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_b,
           bdinfo_file = bdinf_b, locations = locs[6:10])

  bdose_out <- tempfile(fileext = ".bdose"); bdinf_out <- tempfile(fileext = ".bdinfo")
  expect_error(mergebd(bdose_files  = c(bdose_a, bdose_b),
                       bdinfo_files = c(bdinf_a, bdinf_b),
                       bdose_file   = bdose_out,
                       bdinfo_file  = bdinf_out), NA)

  merged <- getbd5info(bdose_out, bdinf_out)

  # All 60 subjects, all 10 SNPs
  expect_equal(nrow(merged$samples), 60L)
  expect_equal(merged$samples$sid,   all_info$samples$sid)
  expect_equal(nrow(merged$snps),    nrow(all_info$snps))
  expect_equal(merged$snps$location, all_info$snps$location)

  # Verify dosage data for all SNPs
  for (i in seq_len(nrow(merged$snps))) {
    orig_snp   <- getsnp(all_info, i, dosageonly = TRUE)
    merged_snp <- getsnp(merged,   i, dosageonly = TRUE)
    expect_equal(merged_snp$dosage, orig_snp$dosage, tolerance = 5e-5,
                 label = paste0("SNP merge dosage SNP ", i))
  }
})

test_that("mergebd SNP merge: subject intersection", {
  bdfile   <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  f5       <- make_bd5(bdfile)
  all_info <- f5$info
  locs     <- all_info$snps$location
  all_sids <- all_info$samples$sid

  # File A: subjects 1-40, SNPs 1-5
  # File B: subjects 11-60, SNPs 6-10
  # Common subjects: 11-40 (30 subjects)
  bdose_a <- tempfile(fileext = ".bdose"); bdinf_a <- tempfile(fileext = ".bdinfo")
  bdose_b <- tempfile(fileext = ".bdose"); bdinf_b <- tempfile(fileext = ".bdinfo")
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_a,
           bdinfo_file = bdinf_a, locations = locs[1:5],
           subjectids = all_sids[1:40])
  subsetbd(bdfiles = c(f5$bdose, f5$bdinf), bdose_file = bdose_b,
           bdinfo_file = bdinf_b, locations = locs[6:10],
           subjectids = all_sids[11:60])

  bdose_out <- tempfile(fileext = ".bdose"); bdinf_out <- tempfile(fileext = ".bdinfo")
  expect_error(mergebd(bdose_files  = c(bdose_a, bdose_b),
                       bdinfo_files = c(bdinf_a, bdinf_b),
                       bdose_file   = bdose_out,
                       bdinfo_file  = bdinf_out), NA)

  merged <- getbd5info(bdose_out, bdinf_out)
  expect_equal(nrow(merged$samples), 30L)
  expect_equal(merged$samples$sid,   all_sids[11:40])
  expect_equal(nrow(merged$snps),    10L)

  # SNPs 1-5: check dosage against original for common subjects
  for (i in 1:5) {
    orig_snp   <- getsnp(all_info, i, dosageonly = TRUE)
    merged_snp <- getsnp(merged,   i, dosageonly = TRUE)
    expect_equal(merged_snp$dosage, orig_snp$dosage[11:40], tolerance = 5e-5,
                 label = paste0("SNP merge intersection dosage SNP ", i))
  }
})
