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
const int HEADERSIZE[NUMFORMATS] = {
  8, 8, 8, 8, 12, 12, 52, 52, 40, 40, 40, 40
};
const int MAXHEADERSIZE = 52;


// Writes the base header for a binary dosage file
// [[Rcpp::export]]
int WriteBinaryDosageBaseHeader(std::string &filename, int formatIndex) {
  std::fstream outfile;
  const int magicWord = 0x65736f62;

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

  outfile.write((char *)&magicWord, sizeof(int));
  outfile.write((char *)(BDFORMATS + formatIndex), sizeof(int));

  outfile.close();
  return 0;
}

//***************************************************************************//
//                        Writing the data                                   //
//***************************************************************************//
// These functions appends data to the end of the binary dosage file
// NOTE: There is no error checking done here. It is assumed this was done
// prior to calling these routines

//  ************************ Constants **************************************//

// Values that dosages and genetic probabilities are multiplied by to change
// them to short integers

const int NUMBEROFBASES = 3;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfff3 is 65,534 or 2^32 - 1
// 0x2710 is 10,000
const unsigned short USBASE[NUMBEROFBASES] = {
  0x7ffe, // Used for format 1.1
  0xfffe, // Used for format 1.2
  0x2710  // Used for all other formats
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
const double DBASE[NUMBEROFBASES] = {
  1. / USBASE[0],
  1. / USBASE[1],
  1. / USBASE[2]
};

// Routine to convert double values to short
// Parameter x - is the vector of doubles to convert
// Parameter us - vector of unsigned shorts to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff
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
      // Missing
      *ps1 = 0xffff;
    } else {
      r1 = (unsigned short)(x[i] * USBASE[base]);
      // The following section checks if r1 or (r1 -1) or (r1 + 1)
      // gives the closest value to the double passed
      x1 = r1 * DBASE[base];
      if (r1 * DBASE[base] < x[i])
        r2 = r1 + 1;
      else
        r2 = r1 - 1;
      x2 = r2 * DBASE[base];

      *ps1 = (fabs(x[i] - x1) < fabs(x[i] - x2)) ? r1 : r2;
    }
    if (i < 10)
      Rcpp::Rcout << x[1] << '\t'
                  << *ps1 << '\t'
                  << r1 << '\t'
                  << x1 << '\t'
                  << r2 << '\t'
                  << x2 << std::endl;
  }
}

// Write only the dosage data
// Paramater filename - name of bindary dosage file
// Parameter dosage - vector of dosage values to write
// Parameter usdosage - vector used to store the converted dosage values. This
//                      is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
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

// Write only the p1 and p2 data
// Paramater filename - name of bindary dosage file
// Parameter p1 - vector of P(g=1) to write
// Parameter p2 - vector of P(g=2) to write
// Parameter usp1 - vector used to store the converted p1 values. This
//                  is passed to avoid allocating and dealllocating memory
// Parameter usp2 - vector used to store the converted p1 values. This
//                  is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
// [[Rcpp::export]]
int WriteBinaryP1P2Data(const std::string &filename,
                            Rcpp::NumericVector &p1,
                            Rcpp::NumericVector &p2,
                            Rcpp::IntegerVector &usp1,
                            Rcpp::IntegerVector &usp2,
                            int base) {
  std::ofstream outfile;

  // Opens file and truncates to size 0. Should already be of size 0.
  outfile.open(filename.c_str(),
               std::ios_base::out | std::ios_base::binary | std::ios_base::app);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  DoubleToUShort(p1, usp1, base);
  DoubleToUShort(p2, usp2, base);

  outfile.close();
  return 0;
}
