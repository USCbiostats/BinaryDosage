#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

const std::string VCFSubjectLineHeader =
  "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

// [[Rcpp::export]]
Rcpp::List OpenVCFFile(Rcpp::StringVector &rFilename) {
  std::string cFilename;
  std::ifstream infile;
  Rcpp::List retVal;
  int subjectLine;
  std::string readBuffer;
  std::streampos dataBegin;
  std::istringstream iss;
  int numSub;
  std::string subId;
  std::vector<std::string> sid, fid;
  Rcpp::DataFrame samples;
  int numSNPs;
  std::streampos lastPos, currentPos;
  std::vector<int> indices;
  int i;

  cFilename = rFilename[0];
  infile.open(cFilename.c_str());
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open VCF file" << std::endl;
    return retVal;
  }
  subjectLine = 0;
  while(1) {
    ++subjectLine;
    std::getline(infile, readBuffer);
    if (infile.fail()) {
      Rcpp::Rcerr << "Unable to find subject line - eof" << std::endl;
      return retVal;
    }
    if (readBuffer[0] != '#') {
      Rcpp::Rcerr << "Unable to find subject line - format error" << std::endl;
      return retVal;
    }
    if (readBuffer[1] != '#') {
      if (readBuffer.substr(1, VCFSubjectLineHeader.size()) != VCFSubjectLineHeader) {
        Rcpp::Rcerr << "Unable to find subject line - column names" << std::endl;
        return retVal;
      }
      break;
    }
  }
  dataBegin = infile.tellg();

  iss.str(readBuffer.substr(VCFSubjectLineHeader.size() + 1));
  numSub = 0;
  while(1) {
    iss >> subId;
    if (iss.fail())
      break;
    ++numSub;
  }
  sid.resize(numSub);
  fid.resize(numSub);
  iss.clear();
  iss.str(readBuffer.substr(VCFSubjectLineHeader.size() + 1));
  for (i = 0; i < numSub; ++i)
    iss >> sid[i];

  numSNPs = 0;
  while(1) {
    std::getline(infile, readBuffer);
    if (infile.fail())
      break;
    ++numSNPs;
  }
  indices.resize(numSNPs);
  infile.clear();
  infile.seekg(dataBegin);
  lastPos = 0;
  for (i = 0; i < numSNPs; ++i) {
    currentPos = infile.tellg();
    indices[i] = (int)(currentPos - lastPos);
    lastPos = currentPos;
    std::getline(infile, readBuffer);
  }

  infile.close();

  samples = Rcpp::DataFrame::create(Rcpp::Named("FID") = fid,
                                    Rcpp::Named("SID") = sid);
  retVal = Rcpp::List::create(Rcpp::Named("filename") = rFilename,
                              Rcpp::Named("headersize") = subjectLine,
                              Rcpp::Named("NumSamples") = numSub,
                              Rcpp::Named("Samples") = samples,
                              Rcpp::Named("Indices") = indices);
  retVal.attr("class") = "genetic-file-info";

  return retVal;
}
