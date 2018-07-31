// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <iostream>
#include <fstream>
#include <Rcpp.h>
#include "GeneticData.h"
#include "BinaryDosage.h"

void WriteMergeHeader(std::ofstream &outfile, int numSubjects, int numSNPs, int numGroups) {
  char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x2 };
  const int zeroInt = 0x0;

  outfile.write(header, 8);
  outfile.write((char *)&numSubjects, sizeof(int));
  outfile.write((char *)&numSNPs, sizeof(int));
  outfile.write((char *)&numGroups, sizeof(int));
  outfile.write((char *)&zeroInt, sizeof(int));
  outfile.write((char *)&zeroInt, sizeof(int));
  outfile.write((char *)&zeroInt, sizeof(int));
  outfile.write((char *)&zeroInt, sizeof(int));
  outfile.write((char *)&zeroInt, sizeof(int));
  outfile.write((char *)&numSubjects, sizeof(int));
}

void WriteMergeSubjects(std::ofstream &outfile, const Rcpp::List &bdList) {
  unsigned int ui, uj;
  std::vector<std::string> subID;

//  for (ui = 0; ui < bdList.size(); ++ui) {
//    subID = Rcpp::as<std::vector<std::string> >(Rcpp::as<Rcpp::List>(bdList[ui])[""])
//    for (uj = 0; uj < Rcpp::as<Rcpp::List>(bdList[ui]).size(); ++uj)
//      outfile.write()
//  }

}

CGeneticData *OpenBinaryDosageFile(const Rcpp::List &geneticInfo) {
  CGeneticData *geneticData = NULL;
  int format, subversion, numSubjects, numSNPs;
  std::string filename;

  format = (int)geneticInfo["format"];
  subversion = (int)geneticInfo["version"];
  numSubjects = (int)geneticInfo["numSubjects"];
  numSNPs = (int)geneticInfo["numSNPs"];
  filename = Rcpp::as<std::string>(geneticInfo["filename"]);

  if (format == 1 && subversion == 1)
    geneticData = new CBinaryDosageFormat1_1(filename, numSubjects, numSNPs);
  else if (format == 1 && subversion == 2)
    geneticData = new CBinaryDosageFormat1_2(filename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 1)
    geneticData = new CBinaryDosageFormat2_1(filename, numSubjects, numSNPs);
  else if (format == 2 && subversion == 2)
    geneticData = new CBinaryDosageFormat2_2(filename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 1)
    geneticData = new CBinaryDosageFormat3_1(filename, numSubjects, numSNPs);
  else if (format == 3 && subversion == 2)
    geneticData = new CBinaryDosageFormat3_2(filename, numSubjects, numSNPs);
  else if (format == 4 && subversion == 1)
    geneticData = new CBinaryDosageFormat4_1(filename, numSubjects, numSNPs);
  else if (format == 4 && subversion == 2)
    geneticData = new CBinaryDosageFormat4_2(filename, numSubjects, numSNPs);
  else
    return NULL;

  return geneticData;
}

//' Function to test opening an array of binary dosage files
//'
//' Function to test opening an array of binary dosage files
//'
//' @param bdList
//' List of binary dosage info lists
//' @param mergeFilename
//' Name of file to merge data into
//' @param snpLoc
//' Locations of snps in the various binary dosage files
//' @return
//' 0 Success
//' 1 failure
//' @export
// [[Rcpp::export]]
int OpenBDFiles(const Rcpp::List &bdList, const std::string &mergeFilename, const Rcpp::NumericMatrix &snpLoc) {
  std::vector<CGeneticData *> bdFiles;
  std::vector<double> d, p0, p1, p2;
  int numSubjects;
  std::ofstream outfile;
  int numFiles;
  int i;

  numFiles = bdList.size();
  bdFiles.resize(numFiles);

  numSubjects = 0;
  for (i = 0; i < numFiles; ++i) {
    bdFiles[i] = OpenBinaryDosageFile(Rcpp::as<Rcpp::List>(bdList[i]));
    numSubjects += bdFiles[i]->NumSubjects();
  }

  outfile.open(mergeFilename.c_str(), std::ios_base::out | std::ios_base::binary);
  Rcpp::Rcout << numSubjects << '\t' << snpLoc.nrow() << std::endl;
  WriteMergeHeader(outfile, numSubjects, snpLoc.nrow(), 1);

  outfile.close();
  for (i = 0; i < numFiles; ++i) {
    if (bdFiles[i])
      delete bdFiles[i];
  }

  return 0;
}

//' Function to test getting dosage values
//'
//' Function to test getting dosage values
//'
//' @param bdInfo
//' Binary dosage info list
//' @return
//' List of dosages
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector GetDosages(const Rcpp::List &bdInfo, int snpNumber) {
  CGeneticData *bdFile = NULL;
  Rcpp::NumericVector d((int)bdInfo["numSubjects"]);
  int i;

  bdFile = OpenBinaryDosageFile(bdInfo);
  if (bdFile == NULL)
    return NULL;

  bdFile->GetFirst();
  for (i = 0; i < bdFile->NumSubjects(); ++i)
    d[i] = bdFile->Dosages()[i];

  if (bdFile)
    delete bdFile;

  return d;
}
