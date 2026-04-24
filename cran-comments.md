## Resubmission

This package was archived on 2025-12-25 because the `src/Makevars` file
specified `CXX_STD = CXX11`, which has been a check NOTE since January 2023.
This specification has been removed; the package now builds using the default
C++17.

## Changes in this version

* Removed `CXX_STD = CXX11` from `src/Makevars`
* Added Format 5 binary dosage files with per-SNP gzip compression
  (`vcftobd5()`, `getbd5info()`, `getbd5snp()`)
* Added `vcftobd()` as the new recommended VCF conversion function
  (uses the vcfppR package); renamed the previous `vcftobd()` to
  `vcftobdlegacy()` with a deprecation message
* Added `updatebd()` to convert any legacy format (1–4) file to Format 5
* Added `subsetbd()` to subset any binary dosage file by MAF, SNP location,
  location range, and/or subject IDs; outputs Format 5
* Added `mergebd()` to merge two or more Format 5 files by subjects or SNPs
* Extended `getsnp()` and `bdapply()` to support Format 5 files

## Test environments

* Windows 11, R 4.5.2
* Linux (Ubuntu), R 4.5.x

## R CMD check results

0 errors | 0 warnings | 0 notes

The Linux check environment produced two environmental notes unrelated to the
package:

* `'qpdf' is needed for checks on size reduction of PDFs` — qpdf is not
  installed on the test server; CRAN machines have qpdf installed.
* `Compilation used the following non-portable flag(s): '-march=x86-64-v3'` —
  this flag comes from the user-level `~/.R/Makevars` on the test server, not
  from the package; `src/Makevars` is empty.
