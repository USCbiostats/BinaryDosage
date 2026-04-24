# BinaryDosage 2.0.0

* Removed `CXX_STD = CXX11` from `src/Makevars`; the package now builds with
  the default C++17 standard.
* Added Format 5 binary dosage files, a new format using per-SNP gzip
  compression and a companion RDS metadata file (`.bdinfo`).
  * `vcftobd5()` — convert a bgzipped VCF file to Format 5.
  * `getbd5info()` — load and validate a Format 5 file pair.
  * `getbd5snp()` — read a single SNP from a Format 5 file by index or ID.
* Added `vcftobd()` as the new recommended function for converting VCF files
  to binary dosage format (requires the vcfppR package). The previous
  `vcftobd()` has been renamed to `vcftobdlegacy()` and emits a deprecation
  message.
* Added `updatebd()` to convert any legacy binary dosage file (formats 1–4)
  to Format 5.
* Added `subsetbd()` to create a new Format 5 file containing a subset of SNPs
  and/or subjects from any binary dosage file (formats 1–5). Filtering options:
  minimum MAF, SNP location vector, SNP location range, and subject ID vector.
* Added `mergebd()` to merge two or more Format 5 files into a single Format 5
  file. Automatically performs a subject merge (unique subjects, common SNPs)
  or a SNP merge (unique SNPs, common subjects) depending on the input files.
* Extended `getsnp()` and `bdapply()` to support Format 5 files.

# BinaryDosage 1.0.0

First release version.

# BinaryDosage 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
