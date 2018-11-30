#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
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
  std::vector<int> index = Rcpp::as< std::vector<int> >(indices);
  std::vector<int>::iterator intIt;
  std::vector<std::streampos> posIndex;
  std::vector<std::streampos>::iterator spIt;
  std::streampos snpIndex;
  int format, version;
  int i, j, n;

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

  snpIndex = 0;
  if (index.size() != 0) {
    posIndex.resize(index.size());
    for(intIt = index.begin(), spIt = posIndex.begin(); intIt != index.end(); ++intIt, ++spIt) {
      snpIndex += *intIt;
      *spIt = snpIndex;
    }
  }

  for (i = 0; i < snpVec.length(); ++i) {
    n = snpVec[i];
    if (n < 1) {
      Rcpp::Rcerr << "SNP number must be a postive value" << std::endl;
      delete miniReader;
      return 1;
    }
    if (n > numSNPs) {
      Rcpp::Rcerr << "SNP index greater than number of SNPs" << std::endl;
      delete miniReader;
      return 1;
    }

    if (snpIndex > 0) {
      if (!miniReader->GetSNP(n, posIndex[n - 1])) {
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

// [[Rcpp::export]]
int GetAlleleFreqC(const std::string &filename, const std::string &filetype,
                   const Rcpp::IntegerVector &subVec, const Rcpp::IntegerVector snpVec,
                   const Rcpp::IntegerVector &indices, Rcpp::NumericVector &freqVec,
                   const int numSubjects, const int numSNPs, const int batchSize) {
  CMiniReader *miniReader = NULL;
  std::vector<int> index = Rcpp::as< std::vector<int> >(indices);
  std::vector<int>::iterator intIt;
  std::vector<std::streampos> posIndex;
  std::vector<std::streampos>::iterator spIt;
  std::streampos snpIndex;
  int format, version;
  std::vector<std::vector<double> > dosages;
  std::vector<std::vector<double> >::iterator ddIt, ddIt2;
  double *d;
  bool lastGroup;
  int firstSNP, lastSNP;
  int i, j, n;

  if (filetype == "BinaryDosage") {
    if (GetBDoseFormat(filename, format, version)) {
      Rcpp::Rcerr << "Unable to open file" << std::endl;
      return 1;
    }
    if (format == 4)
      miniReader = new CBDoseMiniReader4(filename);
    else
      miniReader = new CBDoseMiniReader1(filename, numSubjects, numSNPs);
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

  snpIndex = 0;
  if (index.size() != 0) {
    posIndex.resize(index.size());
    for(intIt = index.begin(), spIt = posIndex.begin(); intIt != index.end(); ++intIt, ++spIt) {
      snpIndex += *intIt;
      *spIt = snpIndex;
    }
  }

  dosages.resize(batchSize);
  for (ddIt = dosages.begin(); ddIt != dosages.end(); ++ddIt)
    ddIt->resize(subVec.size());

  lastGroup = false;
  firstSNP = 0;
  do {
    miniReader->OpenFile();
    lastSNP = firstSNP + batchSize;
    if (lastSNP >= snpVec.length()) {
      lastSNP = snpVec.length();
      lastGroup = true;
    }
    ddIt = dosages.begin();
    for (i = firstSNP; i < lastSNP; ++i, ++ddIt) {
      d = ddIt->data();
      n = snpVec[i];
      if (n < 1) {
        Rcpp::Rcerr << "SNP number must be a postive value" << std::endl;
        delete miniReader;
        return 1;
      }
      if (n > numSNPs) {
        Rcpp::Rcerr << "SNP index greater than number of SNPs" << std::endl;
        delete miniReader;
        return 1;
      }

      if (snpIndex > 0) {
        if (!miniReader->GetSNP(n, posIndex[n - 1])) {
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
      for (j = 0; j < subVec.length(); ++j, ++d)
        *d = miniReader->Dosage()[subVec[j] - 1];
    }
    miniReader->CloseFile();
    for (ddIt2 = dosages.begin(), i = firstSNP; ddIt2 != ddIt; ++ddIt2, ++i)
      freqVec[i] = std::accumulate(ddIt2->begin(), ddIt2->end(), 0.) / (2. * subVec.length());
    firstSNP = lastSNP;
  } while (!lastGroup);

  delete miniReader;

  return 0;
}
