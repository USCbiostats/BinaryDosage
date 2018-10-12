#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "BDoseMiniReader.h"
#include "VCFMiniReader.h"
#include <Rcpp.h>

// [[Rcpp::export]]
int GetSNPValuesC(const std::string &filename, const std::string &filetype, int geneProb,
                  const Rcpp::IntegerVector &subVec, const Rcpp::IntegerVector snpVec,
                  const Rcpp::IntegerVector &indices, Rcpp::NumericMatrix &valueMatrix,
                  const int numSubjects, const int numSNPs) {

  CMiniReader *miniReader = NULL;
  std::ifstream infile;
  std::vector<int> index = Rcpp::as< std::vector<int> >(indices);
  std::streampos snpIndex = 0;
  Rcpp::NumericVector d, p0, p1, p2;
  Rcpp::DataFrame retVal;
  int format, version;
  int i, j, n;
  unsigned int un;
  std::string junk;

  if (filetype == "BinaryDosage") {
    if (GetBDoseFormat(filename, format, version)) {
      Rcpp::Rcerr << "Unable to open file" << std::endl;
      return 1;
    }
    if (format == 4)
      miniReader = new CBDoseMiniReader4(filename);
    else
      miniReader = new CBDoseMiniReader1(filename, numSubjects, numSNPs);
    if (miniReader->P0().size() == 0 && geneProb == 1) {
      Rcpp::Rcerr << "Gene probabilities are requested but don't exist in file" << std::endl;
      delete miniReader;
      return 1;
    }
  } else if (filetype == "VCF") {
    if (index.size() == 0) {
      Rcpp::Rcerr << "VCF files require indices" << std::endl;
      return 1;
    }
    miniReader = new CVCFMiniReader(filename, numSubjects, numSNPs, index[0]);
  } else {
    Rcpp::Rcerr << "Unknown file type" << std::endl;
    return 1;
  }

  if (!miniReader->good()) {
    Rcpp::Rcerr << "Error opening file" << std::endl;
    delete miniReader;
    return 1;
  }

  for (i = 0; i < snpVec.length(); ++i) {
    n = snpVec[i];
    if (n < 1) {
      Rcpp::Rcerr << "SNP number must be a postive value" << std::endl;
      delete miniReader;
      return 1;
    }
    un = (unsigned int)n;

    if (index.size() != 0) {
      if (un > index.size()) {
        Rcpp::Rcerr << "SNP number greater than number of SNP indices" << std::endl;
        delete miniReader;
        return 1;
      }
      snpIndex = 0;
      for (unsigned ui = 0; ui < un; ++ui)
        snpIndex += index[ui];
    }

    if (snpIndex > 0) {
      if (!miniReader->GetSNP(n, snpIndex)) {
        Rcpp::Rcerr << "Error reading SNP" << std::endl;
        delete miniReader;
        return 1;
      }
    } else {
      if (!miniReader->GetSNP(n)) {
        Rcpp::Rcerr << "Error reading SNP" << std::endl;
        delete miniReader;
        return 1;
      }
    }
    if (geneProb == 0) {
      for (j = 0; j < subVec.length(); ++j) {
        valueMatrix(j, i) = miniReader->Dosage()[subVec[j] - 1];
      }
    } else {
      for (j = 0; j < subVec.length(); ++j) {
        valueMatrix(j, 4*i) = miniReader->Dosage()[subVec[j] - 1];
        valueMatrix(j, 4*i + 1) = miniReader->P0()[subVec[j] - 1];
        valueMatrix(j, 4*i + 2) = miniReader->P1()[subVec[j] - 1];
        valueMatrix(j, 4*i + 3) = miniReader->P2()[subVec[j] - 1];
      }
    }
  }

  delete miniReader;

  return 0;
}
