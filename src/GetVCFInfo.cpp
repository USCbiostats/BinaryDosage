#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <Rcpp.h>

void ReadVCFSubjects(std::vector<std::string> &_subjects, const std::string &_subjectString) {
  std::istringstream iss;
  std::string readSubject;

  _subjects.reserve(2000);
  iss.str(_subjectString);
  iss >> readSubject;
  while (!iss.fail()) {
    _subjects.push_back(readSubject);
    iss >> readSubject;
  }
}
//' Function to get information about a VCF file
//'
//' Function to get information about a VCF file.
//' This information is used to get dosage and genetic
//' probabilities from the file.
//'
//' @param vcfFilename
//' Name of VCF file
//' @return
//' List of information about the VCF file
//' @export
// [[Rcpp::export]]
Rcpp::List GetVCFInfo_C(std::string &vcfFile) {
  std::ifstream infile;
  std::string fileFormatString;
  const std::string ExpectedDataHeader = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  std::string dataRead;
  std::vector<std::string> subjects;
  std::istringstream iss;
  std::string readChromosome, readID, readRef, readAlt, readQual, readFilter, readInfo, readFormat;
  int readPOS, curOffset;
  std::vector<std::string> chromosome, id, ref, alt, qual, filter, info, format;
  std::vector<int> pos;
  std::vector<int> offset;
  Rcpp::List retVal;
  Rcpp::DataFrame snps;
  std::streampos lastPos, currentPos;

  infile.open(vcfFile.c_str());
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open vcf file " << vcfFile << std::endl;
    return retVal;
  }

  getline(infile, dataRead);
  if (!infile.good()) {
    Rcpp::Rcerr << "Error reading form vcf file" << std::endl;
    infile.close();
    return retVal;
  }
  if (dataRead.find("##fileformat=") != 0) {
    Rcpp::Rcerr << "File does not appear to be a vcf file" << std::endl;
    Rcpp::Rcerr << "First line:\t" << dataRead << std::endl;
    infile.close();
    return retVal;
  }

  do {
    getline(infile, dataRead);
    if (!infile.good()) {
      Rcpp::Rcerr << "Error reading form vcf file" << std::endl;
      infile.close();
      return retVal;
    }
    if (dataRead[0] != '#') {
      Rcpp::Rcerr << "Unable to find data header row" << std::endl;
      infile.close();
      return retVal;
    }
    if (dataRead[1] != '#')
      break;
  } while (1);

  if (dataRead.find(ExpectedDataHeader)) {
    Rcpp::Rcerr << "Data header row not in expected format" << std::endl;
    infile.close();
    return retVal;
  }

  ReadVCFSubjects(subjects, dataRead.substr(ExpectedDataHeader.length(), dataRead.length() - ExpectedDataHeader.length()));
  Rcpp::Rcout << subjects.size() << '\t' << subjects[0] << '\t' << subjects[subjects.size() - 1] << std::endl;

  chromosome.reserve(600000);
  id.reserve(600000);
  ref.reserve(600000);
  alt.reserve(600000);
  qual.reserve(600000);
  filter.reserve(600000);
  info.reserve(600000);
  format.reserve(600000);
  offset.reserve(600000);

  lastPos = 0;
  currentPos = infile.tellg();
  getline(infile, dataRead);
  while(infile.good()) {
    iss.str(dataRead);
    iss >> readChromosome >> readPOS >> readID >> readRef >> readAlt >> readQual >> readFilter >> readInfo >> readFormat;
    if (iss.fail()) {
      Rcpp::Rcerr << "Error reading data" << std::endl;
      infile.close();
      return retVal;
    }
    curOffset = (int)(currentPos - lastPos);
    chromosome.push_back(readChromosome);
    pos.push_back(readPOS);
    id.push_back(readID);
    ref.push_back(readRef);
    alt.push_back(readAlt);
    qual.push_back(readQual);
    filter.push_back(readFilter);
    info.push_back(readInfo);
    format.push_back(readFormat);
    offset.push_back(curOffset);
    lastPos = currentPos;
    currentPos = infile.tellg();
    getline(infile, dataRead);
  }
  Rcpp::Rcout << chromosome.size() << '\t' << chromosome[0] << '\t' << pos[chromosome.size() - 1] << std::endl;

  infile.close();

  snps = Rcpp::DataFrame::create(Rcpp::Named("Chromosome") = chromosome,
                                 Rcpp::Named("Position") = pos,
                                 Rcpp::Named("ID") = id,
                                 Rcpp::Named("Reference") = ref,
                                 Rcpp::Named("Alternate") = alt,
                                 Rcpp::Named("Quality") = qual,
                                 Rcpp::Named("Filter") = filter,
                                 Rcpp::Named("Information") = info,
                                 Rcpp::Named("Format") = format,
                                 Rcpp::Named("Offset") = offset);
  retVal = Rcpp::List::create(Rcpp::Named("Samples") = subjects,
                              Rcpp::Named("SNPs") = snps);
  return retVal;
}
