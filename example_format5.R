# Example: Converting a VCF file to BinaryDosage Format 5
#
# Format 5 stores each SNP as an individual gzip-compressed block in a
# .bdose file.  A companion .bdinfo file (an RDS list of class
# "genetic-info") holds the sample IDs, SNP table, and byte offsets.
#
# Prerequisites:
#   install.packages("vcfppR")
#
# Run from the package root directory:
#   Rscript example_format5.R
#
# Note: format5.R is sourced directly here because the package has not yet
# been rebuilt with the Format 5 functions.  Once the package is installed,
# replace the two lines below with:  library(BinaryDosage)

library(BinaryDosage)   # loads extdata paths and existing helpers
source("R/format5.R")   # loads vcftobd, getbd5info, getbd5snp

# ---------------------------------------------------------------------------
# 1. Locate the example VCF file included with the package
# ---------------------------------------------------------------------------
vcf_file <- system.file("extdata", "set1a.vcf.gz", package = "BinaryDosage")

# Output files (written to a temporary directory so nothing is left behind)
bdose_file  <- file.path(tempdir(), "set1a.bdose")
bdinfo_file <- file.path(tempdir(), "set1a.bdinfo")

# ---------------------------------------------------------------------------
# 2. Convert the VCF to Format 5
# ---------------------------------------------------------------------------
vcftobd(vcffile  = vcf_file,
         bdose_file  = bdose_file,
         bdinfo_file = bdinfo_file)

cat("Output files created:\n")
cat("  .bdose :", bdose_file,  " (", file.info(bdose_file)$size,  "bytes )\n")
cat("  .bdinfo:", bdinfo_file, " (", file.info(bdinfo_file)$size, "bytes )\n\n")

# ---------------------------------------------------------------------------
# 3. Load the file information
# ---------------------------------------------------------------------------
bd5 <- getbd5info(bdose_file = bdose_file, bdinfo_file = bdinfo_file)

cat("File information\n")
cat("  Class          :", class(bd5), "\n")
cat("  .bdose file    :", bd5$filename, "\n")
cat("  Uses family ID :", bd5$usesfid, "\n")
cat("  One chromosome :", bd5$onechr, "\n")
cat("  Samples        :", nrow(bd5$samples), "\n")
cat("  SNPs           :", nrow(bd5$snps), "\n\n")

cat("First 5 samples:\n")
print(head(bd5$samples, 5))
cat("\n")

cat("First 5 SNPs:\n")
print(head(bd5$snps, 5))
cat("\n")

# ---------------------------------------------------------------------------
# 4. Read a SNP by index
# ---------------------------------------------------------------------------
snp1 <- getbd5snp(bd5info = bd5, snp = 1L)

cat("SNP 1 --", bd5$snps$snpid[1], "\n")
result <- data.frame(
  SampleID = bd5$samples$sid,
  Dosage   = round(snp1$dosage, 4),
  P_00     = round(snp1$p0,     4),
  P_01     = round(snp1$p1,     4),
  P_11     = round(snp1$p2,     4)
)
print(head(result, 10))
cat("\n")

# ---------------------------------------------------------------------------
# 5. Read a SNP by ID
# ---------------------------------------------------------------------------
snp_id <- bd5$snps$snpid[5]
snp5   <- getbd5snp(bd5info = bd5, snp = snp_id)

cat("SNP by ID --", snp_id, "\n")
cat("  Dosage values (first 10):", round(snp5$dosage[1:10], 4), "\n\n")

# ---------------------------------------------------------------------------
# 6. Compute alternate allele frequency for every SNP
# ---------------------------------------------------------------------------
n_snps <- nrow(bd5$snps)
aaf <- numeric(n_snps)
for (i in seq_len(n_snps)) {
  d <- getbd5snp(bd5, i)$dosage
  aaf[i] <- mean(d, na.rm = TRUE) / 2
}

cat("Alternate allele frequencies:\n")
aaf_table <- data.frame(snpid = bd5$snps$snpid, aaf = round(aaf, 4))
print(aaf_table)

# ---------------------------------------------------------------------------
# 7. Clean up
# ---------------------------------------------------------------------------
unlink(c(bdose_file, bdinfo_file))
cat("\nTemporary files removed.\n")
