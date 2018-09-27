#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GetBDoseInfo.h"
#include <Rcpp.h>

// Converts a C++ vector of double vectors to a Rcpp::NumericMatrix
Rcpp::NumericMatrix VectorVectorToMatrix(const std::vector<std::vector<double> > &_vToDbl) {
  unsigned int rows, columns;
  unsigned int ui, uj;
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
Rcpp::List GetBinaryDosageInfoC(const std::string &bdFilename, const std::string &famFilename, const std::string &mapFilename, const int index) {
  Rcpp::List retVal;
  Rcpp::DataFrame samples, snps, snpInfo;
  int format, version;
  bool usesFamilyID;
  unsigned int numSamples;
  unsigned int numSNPs;
  Rcpp::NumericMatrix x1, x2, x3, x4;

  CBDoseReader *bdr = NULL;

  if (GetBDoseFormat(bdFilename, format, version))
    return retVal;

  if (format == 4)
    bdr = new CBDoseReader4(bdFilename);
  else
    bdr = new CBDoseReader1(bdFilename, famFilename, mapFilename);

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
