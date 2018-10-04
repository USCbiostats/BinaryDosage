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

// Function to get infromation about a binary dosage file
//
// Function to get infromation about a binary dosage file
// @param bdFilename
// Name of binary dosage file
// @param famFilename
// Name of subject data file
// @param mapFilename
// Name of SNP data file
// @param indices
// Index of SNPs returned from GetBDoseInfo
// @param valueMatrix
// Matrix to write values to. This needs to be created in R
// @return
// List with information about the binary dosage file
// @export
// [[Rcpp::export]]
int GetSNPValuesC(const std::string &bdFilename, const Rcpp::IntegerVector &subVec, const Rcpp::IntegerVector snpVec,
                  const Rcpp::IntegerVector &indices, Rcpp::NumericMatrix &valueMatrix) {

  CMiniReader *miniReader = NULL;
  std::vector<int> index = Rcpp::as< std::vector<int> >(indices);
  std::streampos snpIndex = 0;
  Rcpp::NumericVector d, p0, p1, p2;
  Rcpp::DataFrame retVal;
  int format, version;
  int i, j, n;
  unsigned int un;

  miniReader = new CBDoseMiniReader4(bdFilename);
  if (!miniReader->good()) {
    Rcpp::Rcerr << "Error reading file" << std::endl;
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

    GetBDoseFormat(bdFilename, format, version);
    if (format != 4) {
      Rcpp::Rcerr << "GetSNPValues currently only works with format 4.X or greater" << std::endl;
      delete miniReader;
      return 1;
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
    if (version == 1) {
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
int GetVCFSNPValues(const std::string &vcfFilename, const Rcpp::IntegerVector &subVec, const Rcpp::IntegerVector snpVec,
                  const Rcpp::IntegerVector &indices, Rcpp::NumericMatrix &valueMatrix, int startRow,
                  int numSubjects, int numSNPs) {
  std::vector<int> index = Rcpp::as< std::vector<int> >(indices);
  std::streampos startData;
  std::ifstream infile;
  std::string junk;
  int i;

  if (index.size() == 0) {
    if (startRow < 3) {
      Rcpp::Rcerr << "Invalid starting location of data" << std::endl;
      return 1;
    }
    infile.open(vcfFilename.c_str());
    if (!infile.good()) {
      Rcpp::Rcerr << "Unable to open VCF file" << std::endl;
      return 1;
    }
    for (i = 0; i < startRow; ++i)
      getline(infile, junk);
    if (!infile.good()) {
      Rcpp::Rcerr << "Invalid starting location of data" << std::endl;
      infile.close();
      return 1;
    }
    startData = infile.tellg();
  } else {
    startData = index[0];
  }

  CVCFMiniReader vcfReader(vcfFilename, numSubjects, numSNPs, startData);

  vcfReader.GetFirst();
  Rcpp::Rcout << vcfReader.Dosage()[0] << '\t' << vcfReader.Dosage()[numSubjects - 1] << std::endl;
  Rcpp::Rcout << vcfReader.P0()[0] << '\t' << vcfReader.P0()[numSubjects - 1] << std::endl;
  Rcpp::Rcout << vcfReader.P1()[0] << '\t' << vcfReader.P1()[numSubjects - 1] << std::endl;
  Rcpp::Rcout << vcfReader.P2()[0] << '\t' << vcfReader.P2()[numSubjects - 1] << std::endl;

  return 0;
}
