#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List GetVCFHeaderC(std::string &vcfFile) {
  const std::string columnNames[9] = { "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};

  Rcpp::List retVal;
  std::ifstream infile;
  std::string junk;
  std::string format;
  std::istringstream iss;
  std::string colName;
  std::string subID;
  std::vector<std::string> sID, fID;
  std::vector<std::string>::iterator strIt;
  int headerLine;
  int startData;
  int numSub;
  int i;

  infile.open(vcfFile.c_str());
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open VCF file" << std::endl;
    return retVal;
  }

  std::getline(infile, junk);
  if (junk.substr(0, 16) != "##fileformat=VCF") {
    Rcpp::Rcerr << "File does not appear to be a VCF file" << std::endl;
    return retVal;
  }
  format = junk.substr(16);

  headerLine = 1;
  do {
    std::getline(infile, junk);
    ++headerLine;
  } while (infile.good() && junk.substr(0,2) == "##");
  if (!infile.good()) {
    Rcpp::Rcerr << "Could not find header line" << std::endl;
    infile.close();
    return retVal;
  }

  iss.str(junk.substr(1));
  for (i = 0; i < 9; ++i) {
    iss >> colName;
    if (colName != columnNames[i]) {
      Rcpp::Rcerr << "Error reading column names" << std::endl;
      infile.close();
      return retVal;
    }
  }
  startData = (int)infile.tellg();
  infile.close();

  iss >> subID;
  numSub = 0;
  while(1) {
    ++numSub;
    iss >> subID;
    if(iss.fail())
      break;
  }

  sID.resize(numSub);
  fID.resize(numSub);
  iss.clear();
  iss.seekg(0);
  for (i = 0; i < 9; ++i)
    iss >> colName;
  for (strIt = sID.begin(); strIt != sID.end(); ++strIt) {
    iss >> subID;
    *strIt = subID;
  }

  retVal = Rcpp::List::create(Rcpp::Named("filetype") = "VCF",
                              Rcpp::Named("filename") = vcfFile,
                              Rcpp::Named("format") = format,
                              Rcpp::Named("NumSamples") = numSub,
                              Rcpp::Named("headerLine") = headerLine,
                              Rcpp::Named("StartData") = startData,
                              Rcpp::Named("Samples") = Rcpp::DataFrame::create(Rcpp::Named("FID") = fID,
                                          Rcpp::Named("SID") = sID,
                                          Rcpp::Named("stringsAsFactors") = false));
  return retVal;
}

// [[Rcpp::export]]
Rcpp::List GetVCFSNPInfoC(std::string &filename, int startData, int reserve) {
  Rcpp::List retVal;
  std::ifstream infile;
  std::vector<std::string> chromosome, snpID, ref, alt, qual, filter, info, format;
  std::vector<int> pos, index;
  // values read in - too lazy to type meaningful names
  std::string s1, s2, s3, s4, s5, s6, s7, s8;
  int i1, offset;
  std::streampos curPos, lastPos;
  std::string junk;
  int numSNPs;

  infile.open(filename.c_str());
  if (!infile.good()) {
    Rcpp::Rcerr << "Error opening VCF file" << std::endl;
    return retVal;
  }

  infile.seekg(startData);
  if (!infile.good()) {
    Rcpp::Rcerr << "Error going to starting point" << std::endl;
    infile.close();
    return retVal;
  }

  if (reserve == 0)
    reserve = 1000000;

  chromosome.reserve(reserve);
  snpID.reserve(reserve);
  ref.reserve(reserve);
  alt.reserve(reserve);
  qual.reserve(reserve);
  filter.reserve(reserve);
  info.reserve(reserve);
  format.reserve(reserve);
  pos.reserve(reserve);
  index.reserve(reserve);

  curPos = infile.tellg();
  infile >> s1 >> i1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8;
  if (infile.fail()) {
    Rcpp::Rcerr << "Error reading first line of data" << std::endl;
    infile.close();
    return retVal;
  }
  lastPos = 0;
  numSNPs = 0;
  while(1) {
    ++numSNPs;
    offset = (int)(curPos - lastPos);
    chromosome.push_back(s1);
    pos.push_back(i1);
    snpID.push_back(s2);
    ref.push_back(s3);
    alt.push_back(s4);
    qual.push_back(s5);
    filter.push_back(s6);
    info.push_back(s7);
    format.push_back(s8);
    index.push_back(offset);
    std::getline(infile, junk);
    lastPos = curPos;
    curPos = infile.tellg();
    infile >> s1 >> i1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8;
    if (infile.fail()) {
      if (infile.eof())
        break;
      Rcpp::Rcerr << "Error reading data from file" << std::endl;
      infile.close();
      return retVal;
    }
  }

  retVal = Rcpp::List::create(Rcpp::Named("NumSNPs") = numSNPs,
                              Rcpp::Named("SNPs") = Rcpp::DataFrame::create(Rcpp::Named("SNPID") = snpID,
                                          Rcpp::Named("Chromosome") = chromosome,
                                          Rcpp::Named("Location") = pos,
                                          Rcpp::Named("Reference") = ref,
                                          Rcpp::Named("Alternate") = alt,
                                          Rcpp::Named("stringsAsFactors") = false),
                              Rcpp::Named("SNPInfo") = Rcpp::DataFrame::create(Rcpp::Named("Qual") = qual,
                                          Rcpp::Named("Filter") = filter,
                                          Rcpp::Named("Info") = info,
                                          Rcpp::Named("Format") = format,
                                          Rcpp::Named("stringsAsFactors") = false),
                              Rcpp::Named("Indices") = index);

  infile.close();
  return retVal;
}
