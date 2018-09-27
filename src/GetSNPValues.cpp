#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "BDoseMiniReader.h"
#include <Rcpp.h>

//' Function to get infromation about a binary dosage file
//'
//' Function to get infromation about a binary dosage file
//' @param bdFilename
//' Name of binary dosage file
//' @param famFilename
//' Name of subject data file
//' @param mapFilename
//' Name of SNP data file
//' @return
//' List with information about the binary dosage file
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame GetSNPValues(const std::string &bdFilename, const int n) {
  CBDoseMiniReader *bdmr = NULL;
  Rcpp::NumericVector d, p0, p1, p2;
  Rcpp::DataFrame retVal;
  int format, version;

  GetBDoseFormat(bdFilename, format, version);
  if (format != 4) {
    Rcpp::Rcerr << "GetSNPValues currently only works with format 4.X or greater" << std::endl;
    return retVal;
  }

  bdmr = new CBDoseMiniReader4(bdFilename);
  if (!bdmr->good()) {
    Rcpp::Rcerr << "Error reading file" << std::endl;
    delete bdmr;
    return retVal;
  }
  if (!bdmr->GetSNP(n)) {
    Rcpp::Rcerr << "Error reading SNP" << std::endl;
    delete bdmr;
    return retVal;
  }
  if (version == 1) {
    retVal = Rcpp::DataFrame::create(Rcpp::Named("Dosage") = bdmr->Dosage());
  } else {
    retVal = Rcpp::DataFrame::create(Rcpp::Named("Dosage") = bdmr->Dosage(),
                                Rcpp::Named("P0") = bdmr->P0(),
                                Rcpp::Named("P1") = bdmr->P1(),
                                Rcpp::Named("P2") = bdmr->P2());

  }
  delete bdmr;

  return retVal;
}
