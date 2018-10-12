library(GxEScanR)
library(qqman) # QQ and manhattan plots

cov <- read.table("example_pheno.txt", header = T)
bdose <- GxEScanR::GetBinaryDosageInfo("example.bdose")
GxEScan(cov, bdose, "example_results.out")


# Plot Results
results <- read.table("example_results.out", header = T)

results$zP <- 2*pnorm(-abs(results$zG))
results$zGxEP <- 2*pnorm(-abs(results$zGxE))

qq(results$zP)
manhattan(results, chr= 'CHR', bp = 'BP', p = 'zP')

qq(results$zGxEP)
manhattan(results, chr= 'CHR', bp = 'BP', p = 'zGxEP')



