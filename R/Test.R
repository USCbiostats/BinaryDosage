# Create the binary dosage file - only needs to be run once.
# This takes about 4 minutes on a file with about 2500 subjects and 500,000 SNPs on my desktop
VCF2BD("chr21test.dose.vcf", "chr21test.info", "chr21test.bdose")

# Function that will be evaluated on each SNP. This calculates the alternate
# allele frequency based on the dosage values
# OutputData must be a data.table to allow its values to be changed
AlleleFrequency <- function(DosageData, OutputData, snpsUsed, firstCol, lastCol, otherData = NULL) {
  OutputData[snpsUsed, dosefreq := colMeans(DosageData[,firstCol:lastCol]) / 2]
}

# maximum number of SNPs to read in at one time. Currently 200 is the limit
# I need to learn more about data.tables to increase this limit.
maxSNPs <- 200
# Reads in the information about the subjects and SNPs in the binary dosage file
chr21TestInfo <- GetBinaryFileInfo("chr21test.bdose")
# Select the SNPs that will be passed to the function
# The first line subsets the SNPs based on r-squared > 0.7 and allele
# frequency between 0.05 and 0.95. The second line selects all SNPs.
#### snps2use <- which((chr21TestInfo$SNPs$rSquared > 0.7) & (chr21TestInfo$SNPs$AltFreq > 0.05) & (chr21TestInfo$SNPs$AltFreq < 0.95))
snps2use <- 1:chr21TestInfo$FileData$NumSNPs
# This is the data.table that will be used to store the results
# The column dosefreq is where the alternate allele frequency
# calculated from the dosage values will be saved
# Notice that the AlleleFrequency function above using dosefreq
snpFreq <- data.table(chr21TestInfo$SNPs[c(snps2use),])
snpFreq[, dosefreq := as.numeric(0)]
# Calculate the frequencies using the dosage values for the selected SNPs
SNPApply(binaryDosageInfo = chr21TestInfo, function2apply = AlleleFrequency, outputValues = snpFreq, snps2use = snps2use, maxSNPs = maxSNPs)
# View the results
View(snpFreq)
