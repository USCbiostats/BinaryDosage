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
// 0x7ffe is 32,767 or 2^16 - 1
// 0x2710 is 10,000
const unsigned short USBASE[2] = {
  0x7ffe,
  0x2710
};
const double DBASE[2] = {
  1. / USBASE[0],
  1. / USBASE[1]
};

void DoubleToUShort(Rcpp::NumericVector &x,
                    Rcpp::IntegerVector &us,
                    const int base) {
  unsigned short r1, r2;
  unsigned short *ps1;
  double x1, x2;
  int i;

  r1 = r2 = 0;
  x1 = x2 = 0.;
  ps1 = (unsigned short *)&us[0];
  for (i = 0; i < x.size(); ++i, ++ps1) {
    if (x[i] != x[i]) {
      *ps1 = 0xffff;
    } else {
      r1 = (unsigned short)(x[i] * USBASE[base]);
      x1 = r1 * DBASE[base];
      if (r1 * DBASE[base] < x[i])
        r2 = r1 + 1;
      else
        r2 = r1 - 1;
      x2 = r2 * DBASE[base];

      *ps1 = (fabs(x[i] - x1) < fabs(x[i] - x2)) ? r1 : r2;
    }

    if (i < 10)
      Rcpp::Rcout << x[i] << '\t'
                  << *ps1 << '\t'
                  << r1 << '\t'
                  << x1 << '\t'
                  << r2 << '\t'
                  << x2 << std::endl;
  }
}

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
  std::fstream outfile;
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

// [[Rcpp::export]]
int WriteBinaryDosageData(const std::string &filename,
                          Rcpp::NumericVector &dosage,
                          Rcpp::IntegerVector &usdosage,
                          int base) {
  std::ofstream outfile;

  // Opens file and truncates to size 0. Should already be of size 0.
  outfile.open(filename.c_str(),
               std::ios_base::out | std::ios_base::binary | std::ios_base::app);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  DoubleToUShort(dosage, usdosage, base);

  outfile.close();
  return 0;
}

// [[Rcpp::export]]
int WriteBinaryDosageP1Data(const std::string &filename,
                            Rcpp::NumericVector &dosage,
                            Rcpp::NumericVector &p1,
                            Rcpp::IntegerVector &usdosage,
                            Rcpp::IntegerVector &usp1,
                            int base) {
  std::ofstream outfile;

  // Opens file and truncates to size 0. Should already be of size 0.
  outfile.open(filename.c_str(),
               std::ios_base::out | std::ios_base::binary | std::ios_base::app);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  DoubleToUShort(dosage, usdosage, base);
  DoubleToUShort(p1, usp1, base);

  outfile.close();
  return 0;
}
