#include <Rcpp.h>
#include <iostream>
#include <fstream>

// [[Rcpp::export]]
int WriteBinaryDosageHeader(std::string &filename, int format, int subformat) {
  std::ofstream outfile;
  const char magicWord[4] = { 'b', 'o', 's', 'e'};
  char formatString[4] = {0x0, 0x0, 0x0, 0x0};
  int headerSize;
  const char zeroChar[1] = {0x0};

  outfile.open(filename.c_str());
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }
  outfile.close();

  outfile.open(filename.c_str(),
               std::ios_base::out | std::ios_base::in | std::ios_base::binary | std::ios_base::trunc);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  switch(format) {
  case 1:
    headerSize = 8;
    formatString[1] = 0x1;
    break;
  case 2:
    headerSize = 8;
    formatString[1] = 0x2;
    break;
  case 3:
    formatString[1] = 0x3;
    if (subformat < 3)
      headerSize = 12;
    else
      headerSize = 40;
    break;
  case 4:
    formatString[1] = 0x4;
    if (subformat < 3)
      headerSize = 40;
    else
      headerSize = 24;
    break;
  default:
    Rcpp::Rcerr << "Unknown format" << std::endl;
    outfile.close();
    return 1;
    break;
  }
  switch(subformat) {
  case 1:
    formatString[3] = 0x1;
    break;
  case 2:
    formatString[3] = 0x2;
    break;
  case 3:
    formatString[3] = 0x3;
    break;
  case 4:
    formatString[3] = 0x4;
    break;
  default:
    Rcpp::Rcerr << "Unknown subformat" << std::endl;
    outfile.close();
    return 1;
    break;
  }

  for (int i = 0; i < headerSize; ++i)
    outfile.write(zeroChar, 1);
  outfile.seekp(0);
  outfile.write(magicWord, 4);
  outfile.write(formatString, 4);
  if (outfile.fail()) {
    Rcpp::Rcerr << "Error writing to ouput file" << std::endl;
    outfile.close();
    return 1;
  }

  outfile.close();
  return 0;
}
