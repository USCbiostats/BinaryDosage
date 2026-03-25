## Resubmission

This package was archived on 2025-12-25 because the `src/Makevars` file
specified `CXX_STD = CXX11`, which has been a check NOTE since January 2023.
This specification has been removed; the package now builds using the default
C++17.

## Changes in this version

* Removed `CXX_STD = CXX11` from `src/Makevars`
* Added Format 5 binary dosage files: faster conversion and smaller file sizes
  using gzip compression (new `vcftobd()` function via the vcfppR package)
* Renamed the previous `vcftobd()` to `vcftobdlegacy()` with a deprecation
  warning directing users to the new `vcftobd()`
* Extended `getsnp()` and `bdapply()` to support Format 5 files

## Test environments

* Windows 11, R 4.5.2

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE is from CRAN incoming feasibility checks:
* This is a resubmission of an archived package
* The development version number will be updated to a release version before
  final submission
* Three URLs in README.md have been updated to their current locations
