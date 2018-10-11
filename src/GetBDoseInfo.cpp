#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GetBDoseInfo.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List GetBDose4Header(std::string &filename) {
  char header[4];
  char format[4];
  int numSubjects;
  int numSNPs;
  int numGroups;
  int subjectOptions;
  int snpOptions;
  int startSubjectData;
  int startSNPData;
  int startDosageData;
  std::ifstream infile;
  Rcpp::List retVal;

  infile.open(filename.c_str(), std::ios::in | std::ios::binary);
  if (!infile.good())
    return retVal;

  infile.read(header, 4);
  infile.read(format, 4);
  infile.read((char *)&numSubjects, sizeof(int));
  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&numGroups, sizeof(int));
  infile.read((char *)&subjectOptions, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  infile.read((char *)&startSubjectData, sizeof(int));
  infile.read((char *)&startSNPData, sizeof(int));
  infile.read((char *)&startDosageData, sizeof(int));

  if (!infile.good())
    return retVal;

  retVal = Rcpp::List::create(Rcpp::Named("Subjects") = numSubjects,
                              Rcpp::Named("SNPs") = numSNPs,
                              Rcpp::Named("Groups") = numGroups,
                              Rcpp::Named("SubOption") = subjectOptions,
                              Rcpp::Named("SNPOption") = snpOptions,
                              Rcpp::Named("StartSub") = startSubjectData,
                              Rcpp::Named("StartSNP") = startSNPData,
                              Rcpp::Named("StartDose") = startDosageData);
  infile.close();
  return retVal;
}
// Converts a C++ vector of double vectors to a Rcpp::NumericMatrix
Rcpp::NumericMatrix VectorVectorToMatrix(const std::vector<std::vector<double> > &_vToDbl) {
  int rows, columns;
  int ui, uj;
  std::vector<std::vector<double> >::const_iterator vDblIt;
  std::vector<double>::const_iterator dblIt;

  rows = _vToDbl.size();
  columns = _vToDbl[0].size();
  //  Rcpp::Rcout << rows << '\t' << columns << std::endl;

  Rcpp::NumericMatrix x(rows, columns);
  ui = 0;
  for (vDblIt = _vToDbl.begin(); vDblIt != _vToDbl.end(); ++vDblIt, ++ui) {
    uj = 0;
    for(dblIt = vDblIt->begin(); dblIt != vDblIt->end(); ++ dblIt, ++uj)
      x(ui, uj) = *dblIt;
  }
  return x;
}

// [[Rcpp::export]]
Rcpp::List GetBDoseFormatC(const std::string &bdFilename) {
  int format, version;
  Rcpp::List retVal;

  if (GetBDoseFormat(bdFilename, format, version))
    return retVal;
  retVal = Rcpp::List::create(Rcpp::Named("Format") = format,
                              Rcpp::Named("Version") = version);
  return retVal;
}

// [[Rcpp::export]]
Rcpp::List GetBinaryDosage1Info(const std::string &bdFilename, const Rcpp::DataFrame &subjects,
                                const Rcpp::DataFrame &snps, const int index) {
  Rcpp::List retVal;
  Rcpp::DataFrame samples, SNPs, snpInfo;
  bool usesFamilyID;
  std::vector<std::string> fid = Rcpp::as<std::vector<std::string> >(subjects["FID"]);
  std::vector<std::string> sid = Rcpp::as<std::vector<std::string> >(subjects["SID"]);
  std::vector<std::string> snpID = Rcpp::as<std::vector<std::string> >(snps["SNPID"]);
  std::vector<std::string> chromosome = Rcpp::as<std::vector<std::string> >(snps["CHR"]);
  std::vector<std::string> refAllele = Rcpp::as<std::vector<std::string> >(snps["REF"]);
  std::vector<std::string> altAllele = Rcpp::as<std::vector<std::string> >(snps["ALT"]);
  std::vector<int> location = Rcpp::as<std::vector<int> >(snps["LOC"]);
  CBDoseReader1 bdr(bdFilename, fid, sid, snpID, chromosome, location, refAllele, altAllele);

  if (!bdr.good()) {
    Rcpp::Rcerr << "Unable to open binary dosage file" << std::endl;
    return retVal;
  }

  samples = Rcpp::DataFrame::create(Rcpp::Named("FID") = bdr.FamilyID(),
                                    Rcpp::Named("SID") = bdr.SampleID(),
                                    Rcpp::Named("stringsAsFactors") = false);
  SNPs = Rcpp::DataFrame::create(Rcpp::Named("SNPID") = bdr.SNPID(),
                                 Rcpp::Named("Chromosome") = bdr.Chromosome(),
                                 Rcpp::Named("Location") = bdr.Location(),
                                 Rcpp::Named("Reference") = bdr.ReferenceAllele(),
                                 Rcpp::Named("Alternate") = bdr.AlternateAllele(),
                                 Rcpp::Named("stringsAsFactors") = false);
  snpInfo = Rcpp::DataFrame::create(Rcpp::Named("Dosage") = bdr.Dosage(),
                                    Rcpp::Named("P0") = bdr.P0(),
                                    Rcpp::Named("P1") = bdr.P1(),
                                    Rcpp::Named("p2") = bdr.P2());

  retVal = Rcpp::List::create(Rcpp::Named("filetype") = "BinaryDosage",
                              Rcpp::Named("filename") = bdFilename,
                              Rcpp::Named("format") = bdr.Format(),
                              Rcpp::Named("version") = bdr.Version(),
                              Rcpp::Named("Groups") = bdr.NumGroups(),
                              Rcpp::Named("GroupSizes") = bdr.GroupSize(),
                              Rcpp::Named("NumSamples") = bdr.NumSamples(),
                              Rcpp::Named("usesFID") = usesFamilyID,
                              Rcpp::Named("Samples") = samples,
                              Rcpp::Named("NumSNPs") = bdr.NumSNPs(),
                              Rcpp::Named("SNPs") = SNPs,
                              Rcpp::Named("SNPInfo") = snpInfo,
                              Rcpp::Named("Indices") = bdr.Indices());
  return retVal;
}

// Function to get infromation about a binary dosage file
//
// Function to get infromation about a binary dosage file
// @param bdFilename
// Name of binary dosage file
// @param famFilename
// Name of subject data file
// @param mapFilename
// Name of SNP data file
// @return
// List with information about the binary dosage file
// @export
// [[Rcpp::export]]
Rcpp::List GetBinaryDosage4Info(const std::string &bdFilename, const int index) {
  Rcpp::List retVal;
  Rcpp::DataFrame samples, snps, snpInfo;
  int format, version;
  bool usesFamilyID;
  int numSamples;
  int numSNPs;
  Rcpp::NumericMatrix x1, x2, x3, x4;
  std::ifstream infile;
  CBDoseReader *bdr = NULL;

  if (GetBDoseFormat(bdFilename, format, version))
    return retVal;

  if (format != 4)
    return retVal;

  bdr = new CBDoseReader4(bdFilename);
  if (!bdr->good())
    return retVal;
  if (index)
    bdr->GetIndices();

  if (bdr->FamilyID()[0] == "")
    usesFamilyID = false;
  else
    usesFamilyID = true;
  numSamples = bdr->SampleID().size();
  samples = Rcpp::DataFrame::create(Rcpp::Named("FID") = bdr->FamilyID(),
                                    Rcpp::Named("SID") = bdr->SampleID(),
                                    Rcpp::Named("stringsAsFactors") = false);

  numSNPs = bdr->Location().size();
  snps = Rcpp::DataFrame::create(Rcpp::Named("SNPID") = bdr->SNPID(),
                                 Rcpp::Named("Chromosome") = bdr->Chromosome(),
                                 Rcpp::Named("Location") = bdr->Location(),
                                 Rcpp::Named("Reference") = bdr->ReferenceAllele(),
                                 Rcpp::Named("Alternate") = bdr->AlternateAllele(),
                                 Rcpp::Named("stringsAsFactors") = false);

  x1 = VectorVectorToMatrix(bdr->AlternateAlleleFreq());
  x2 = VectorVectorToMatrix(bdr->MinorAlleleFreq());
  x3 = VectorVectorToMatrix(bdr->AvgCall());
  x4 = VectorVectorToMatrix(bdr->RSquared());

  snpInfo = Rcpp::DataFrame::create(Rcpp::Named("AAF") = x1,
                                    Rcpp::Named("MAF") = x2,
                                    Rcpp::Named("AvgCall") = x3,
                                    Rcpp::Named("RSquared") = x4,
                                    Rcpp::Named("stringsAsFactors") = false);
  retVal = Rcpp::List::create(Rcpp::Named("filetype") = "BinaryDosage",
                              Rcpp::Named("filename") = bdFilename,
                              Rcpp::Named("format") = format,
                              Rcpp::Named("version") = version,
                              Rcpp::Named("Groups") = bdr->NumGroups(),
                              Rcpp::Named("GroupSizes") = bdr->GroupSize(),
                              Rcpp::Named("NumSamples") = numSamples,
                              Rcpp::Named("usesFID") = usesFamilyID,
                              Rcpp::Named("Samples") = samples,
                              Rcpp::Named("NumSNPs") = numSNPs,
                              Rcpp::Named("SNPs") = snps,
                              Rcpp::Named("SNPInfo") = snpInfo,
                              Rcpp::Named("Indices") = bdr->Indices());

  if (bdr)
    delete bdr;

  return retVal;
}
