#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "ReadBinaryDosage.h"
#include <Rcpp.h>

int GetVersion(std::ifstream &infile, std::string &status, int &version, int &subVersion) {
  char header[8];
  const char first4[4] = { 'b', 'o', 's', 'e' };

  status = "Good";
  version = 0;
  subVersion = 0;

  infile.read(header, 8);
  if (!infile.good()) {
    status = "Error reading header";
    return 1;
  }
  if (memcmp(header, first4, 4)) {
    status = "Not a binary dosage file";
    return 1;
  }
  if (header[4] != 0 || header[6] != 0) {
    status = "Unknown format";
    return 1;
  }
  version = header[5];
  subVersion = header[7];
  if (version < 1 || version > 4 || subVersion < 1 || subVersion > 2) {
    status = "Unknown format";
    return 1;
  }
  return 0;
}
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
Rcpp::List GetBinaryDosageInfoC(const std::string &bdFilename, const std::string &famFilename, const std::string &mapFilename) {
  std::string status = "Good";
  int version, subVersion;
  std::ifstream infile;
  CReadBinaryDosageX *bdFile = NULL;
  std::vector<std::string> filenames;

  Rcpp::List bdInfo;
  Rcpp::DataFrame subjects;
  int usesFID = 1;
  std::vector<std::string> famID;
  Rcpp::DataFrame snps;
  std::vector<std::string> snpID, chromosome, refAllele, altAllele;
  int listSize;

  std::vector<std::string>::iterator strIt1, strIt2;
  std::vector<int>::const_iterator intIt;

  infile.open(bdFilename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!infile.good()) {
    status = "Unable to open binary dosage file";
    bdInfo = Rcpp::List::create(Rcpp::Named("Status") = status);
    return bdInfo;
  }
  if (GetVersion(infile, status, version, subVersion)) {
    bdInfo = Rcpp::List::create(Rcpp::Named("Status") = status);
    return bdInfo;
  }
  infile.close();

  if (version < 4) {
    if (famFilename == "" || mapFilename == "") {
      status = "Binary dosage files in formats 1, 2, and 3 require a family and map file";
      bdInfo = Rcpp::List::create(Rcpp::Named("Status") = status);
      return bdInfo;
    }
    filenames.resize(3);
    filenames[0] = bdFilename;
    filenames[1] = famFilename;
    filenames[2] = mapFilename;
  } else {
    if (famFilename != "" || mapFilename != "") {
      status = "Binary dosage files in format 4 or greater do not use a family and map file";
      bdInfo = Rcpp::List::create(Rcpp::Named("Status") = status);
      return bdInfo;
    }
    filenames.resize(1);
    filenames[0] = bdFilename;
  }

  switch (subVersion) {
  case 1:
    switch (version) {
    case 1:
      bdFile = new CReadBinaryDosage11(filenames);
      break;
    case 2:
      bdFile = new CReadBinaryDosage21(filenames);
      break;
    case 3:
      bdFile = new CReadBinaryDosage31(filenames);
      break;
    case 4:
      bdFile = new CReadBinaryDosage41(filenames);
      break;
    default:
      status = "Can't get here";
      break;
    }
    break;
  case 2:
    switch (version) {
    case 1:
      bdFile = new CReadBinaryDosage12(filenames);
      break;
    case 2:
      bdFile = new CReadBinaryDosage22(filenames);
      break;
    case 3:
      bdFile = new CReadBinaryDosage32(filenames);
      break;
    case 4:
      bdFile = new CReadBinaryDosage42(filenames);
      break;
    default:
      status = "Can't get here";
      break;
    }
    break;
  default:
    status = "Can't get here";
    break;
  }

  bdFile->ReadHeader();
  bdFile->ReadSubjects();
  bdFile->ReadGroups();
  bdFile->ReadSNPs();

  if (bdFile->FamilyID().size() == 0) {
    usesFID = 0;
    famID.resize(bdFile->NumSubjects());
    for (strIt1 = famID.begin(); strIt1 != famID.end(); ++strIt1)
      *strIt1 = "";
  } else {
    famID = bdFile->FamilyID();
  }
  subjects = Rcpp::DataFrame::create(Rcpp::Named("FID") = famID,
                                     Rcpp::Named("IID") = bdFile->SubjectID(),
                                     Rcpp::Named("stringsAsFactors") = false);

  if(bdFile->Chromosome().size() == 1) {
    chromosome.resize(bdFile->NumSNPs());
    for (strIt1 = chromosome.begin(); strIt1 != chromosome.end(); ++strIt1)
      *strIt1 = bdFile->Chromosome()[0];
  } else {
    chromosome = bdFile->Chromosome();
  }

  if(bdFile->SNPID().size() == 0) {
    snpID.resize(bdFile->NumSNPs());
    intIt = bdFile->Location().begin();
    for (strIt1 = snpID.begin(), strIt2 = chromosome.begin(); strIt1 != snpID.end(); ++strIt1, ++strIt2, ++intIt)
      *strIt1 = *strIt2 + std::to_string(*intIt);
  } else {
    snpID = bdFile->SNPID();
  }

  if (bdFile->ReferenceAllele().size() == 0) {
    refAllele.resize(bdFile->NumSNPs());
    for (strIt1 = refAllele.begin(); strIt1 != refAllele.end(); ++strIt1)
      *strIt1 = "1";
  } else {
    refAllele = bdFile->ReferenceAllele();
  }

  if (bdFile->AlternateAllele().size() == 0) {
    altAllele.resize(bdFile->NumSNPs());
    for (strIt1 = altAllele.begin(); strIt1 != altAllele.end(); ++strIt1)
      *strIt1 = "2";
  } else {
    altAllele = bdFile->AlternateAllele();
  }

  snps = Rcpp::DataFrame::create(Rcpp::Named("SNP") = snpID,
                                 Rcpp::Named("CHR") = chromosome,
                                 Rcpp::Named("BP") = bdFile->Location(),
                                 Rcpp::Named("A1") = refAllele,
                                 Rcpp::Named("A2") = altAllele,
                                 Rcpp::Named("stringsAsFactors") = false);


  bdInfo = Rcpp::List::create(Rcpp::Named("Status") = status,
                              Rcpp::Named("Format") = version,
                              Rcpp::Named("Version") = subVersion,
                              Rcpp::Named("NumSub") = bdFile->NumSubjects(),
                              Rcpp::Named("UsesFID") = usesFID,
                              Rcpp::Named("Subjects") = subjects,
                              Rcpp::Named("NumSNPs") = bdFile->NumSNPs(),
                              Rcpp::Named("SNPs") = snps);

  if (bdFile)
    delete bdFile;
  return bdInfo;
}
