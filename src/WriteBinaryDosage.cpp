#include <Rcpp.h>
#include <iostream>
#include <fstream>

const int NUMFORMATS = 12;
const int BDFORMATS[NUMFORMATS] = {
  0x01000100,
  0x02000100,
  0x01000200,
  0x02000200,
  0x01000300,
  0x02000300,
  0x03000300,
  0x04000300,
  0x01000400,
  0x02000400,
  0x03000400,
  0x04000400
};
// [[Rcpp::export]]
int BinaryDosageHeaderFormat(int format, int subformat) {
  int bdformat = 0;
  unsigned char &cformat = *((unsigned char *)&bdformat + 1);
  unsigned char &csubformat = *((unsigned char *)&bdformat + 3);

  cformat = (unsigned char)format;
  csubformat = (unsigned char)subformat;
  for (int i = 0; i < NUMFORMATS; ++i) {
    if (BDFORMATS[i] == bdformat)
      return bdformat;
  }

  return 0;
}

// [[Rcpp::export]]
int WriteBinaryDosageHeader(std::string &filename, int format, int subformat) {
  std::ofstream outfile;
  int header[2] = {0x65736f62, 0};
  int &fileformat = header[1];

// Create the file - if file already exists, truncates to size 0.
  outfile.open(filename.c_str());
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }
  outfile.close();

  // Opens file and truncates to size 0. Should already be of size 0.
  outfile.open(filename.c_str(),
               std::ios_base::out | std::ios_base::in | std::ios_base::binary | std::ios_base::trunc);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  fileformat = BinaryDosageHeaderFormat(format, subformat);
  Rcpp::Rcerr << std::hex << fileformat << std::endl;
  outfile.close();
  return 0;
}
