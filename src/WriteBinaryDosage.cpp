#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

extern const std::ios_base::openmode READWRITEBINARY = std::ios_base::in | std::ios_base::out | std::ios_base::binary | std::ios_base::ate;
extern const std::ios_base::openmode READBINARY = std::ios_base::in | std::ios_base::binary;
extern const std::ios_base::openmode WRITEBINARY = std::ios_base::out | std::ios_base::binary | std::ios_base::app;

//***************************************************************************//
//                        Writing the header                                 //
//***************************************************************************//
// These functions write the headers to the binary dosage files
// NOTE: There is no error checking done here. It is assumed this was done
// prior to calling these routines

//  ************************ Constants **************************************//
// Magic word for binary dosage files
extern const int MAGICWORD = 0x65736f62;
// Format ID stored in file
extern const std::vector<std::vector<int> > FORMAT = {
  { 0x01000100, 0x01000200},
  { 0x02000100, 0x02000200},
  { 0x03000100, 0x03000200, 0x03000300, 0x03000400},
  { 0x04000100, 0x04000200, 0x04000300, 0x04000400}
};
// Size of the header for each format in bytes
const std::vector<std::vector<int> > HEADERSIZE = {
  {8, 8},
  {8, 8},
  {12, 12, 72, 72},
  {40, 40, 24, 24}
};

//  ************** Write the base header for all formats*********************//

// Writes the base header for a binary dosage file
// Parameter filename - Name of binary dosage file
// Parameter format - Foramt of the binary dosage file
// Parameter subformat - Subformat of the binary dosage file
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageBaseHeader(std::string &filename, int format, int subformat) {
  std::ofstream outfile;

  // Open the file - if file already exists, truncates to size 0.
  // Only opens for output
  outfile.open(filename.c_str(), std::ios_base::out | std::ios_base::binary);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.write((char *)&MAGICWORD, sizeof(int));
  outfile.write((char *)&FORMAT[format][subformat], sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write the additional information for headers 3.1 and 3.2 *****//

// Writes the additional header info for formats 3.1 and 3.2
// Parameter filename - Name of binary dosage file
// Parameter numSubjects - number of subjects in data
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader3A(std::string &filename, int numSubjects) {
  std::ofstream outfile;

  // Open the file for appending
  // Only opens for output
  outfile.open(filename.c_str(), WRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.write((char *)&numSubjects, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write the additional information for headers 3.3 and 3.4 *****//

// Parameter filename - Name of binary dosage file
// Parameter md5samples - MD5 hash for the samples data frame
// Parameter md5SNPs - MD5 hash for the SNPs data frame
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader3B(std::string &filename,
                              std::string &md5samples,
                              std::string &md5SNPs) {
  std::ofstream outfile;

  // Open the file for appending
  // Only opens for output
  outfile.open(filename.c_str(), WRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.write(md5samples.c_str(), 32);
  outfile.write(md5SNPs.c_str(), 32);

  outfile.close();
  return 0;
}

//  ********** Write the additional information for headers 4.1 and 4.2 *****//

// Writes the additional header info for formats 4.1 and 4.2
// Parameter filename - Name of binary dosage file
// Parameter numSubjects - number of subjects in data
// Parameter numSubjects - number of SNPs in data
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader4A(std::string &filename, int numSubjects, int numSNPs) {
  std::fstream outfile;
  const int zero = 0;

  // Open the file for appending
  // Only opens for output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekp(8);
  // Zero out the rest of the data. It will be filled in later
  for (int i = 0; i < 8; ++i)
    outfile.write((char *)&zero, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write the additional information for headers 4.3 and 4.4 *****//

// Writes the additional header info for formats 4.3 and 4.4
// Parameter filename - Name of binary dosage file
// Return - 0 successful, 1 error
// [[Rcpp::export]]
int WriteBinaryDosageHeader4B (std::string &filename, int numSubjects, int numSNPs) {
  std::ofstream outfile;
  const int zero = 0;

  // Open the file for appending
  // Only opens for output
  outfile.open(filename.c_str(), WRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  // Zero out the data. It will be filled in later
  for (int i = 0; i < 4; ++i)
    outfile.write((char *)&zero, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write group information for formats 4.1 and 4.2 **************//

// [[Rcpp::export]]
int WriteBDGroups(std::string &filename, Rcpp::IntegerVector &groups) {
  std::fstream outfile;
  int numGroups, groupsize;
  int subjectOffset;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekp(16);
  numGroups = groups.length();
  outfile.write((char *)&numGroups, sizeof(int));

  outfile.seekp(40);
  for (int i = 0; i < numGroups; ++i) {
    groupsize = groups[i];
    outfile.write((char *)&groupsize, sizeof(int));
  }

  subjectOffset = (int)outfile.tellp();
  outfile.seekp(28);
  outfile.write((char *)&subjectOffset, sizeof(int));

  outfile.close();

  return 0;
}

//  ********** Write group information for formats 4.3 and 4.4 **************//

// [[Rcpp::export]]
int WriteBDGroups2(std::string &filename, Rcpp::IntegerVector &groups) {
  std::fstream outfile;
  int numGroups, groupsize;
  int subjectOffset;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekp(24);
  numGroups = groups.length();
  outfile.write((char *)&numGroups, sizeof(int));

  for (int i = 0; i < numGroups; ++i) {
    groupsize = groups[i];
    outfile.write((char *)&groupsize, sizeof(int));
  }

  subjectOffset = (int)outfile.tellp();
  outfile.seekp(8);
  outfile.write((char *)&subjectOffset, sizeof(int));

  outfile.close();

  return 0;
}

//  ********** Write a string to the header for format 4 ********************//

int WriteBDString(std::fstream &outfile, std::string &outstring) {
  char zero = 0x0;

  if (outstring.length() > 0) {
    outfile.write(outstring.c_str(), outstring.length());
    outfile.write(&zero, 1);
  }

  return 0;
}

//  ********** Write a integer vector to the header for format 4 ************//

int WriteBDInteger(std::fstream &outfile, Rcpp::IntegerVector &outvector) {

  if (outvector.length() > 0)
    outfile.write((char *)&outvector[0], outvector.length() * sizeof(int));

  return 0;
}

//  ********** Write a numeric vector to the header for format 4 ************//

int WriteBDNumeric(std::fstream &outfile, Rcpp::NumericVector &outvector) {

  if (outvector.length() > 0)
    outfile.write((char *)&outvector[0], outvector.length() * sizeof(double));

  return 0;
}

//  ********** Write a the family information for format 4 ******************//

// [[Rcpp::export]]
int WriteBDFamilyInfoC(std::string &filename,
                       int numSub,
                       std::string &sid,
                       std::string &fid,
                       int numSubLoc,
                       int suboffsetLoc,
                       int snpoffsetLoc) {
  std::fstream outfile;
  int sidsize, fidsize;
  int suboffset, snpoffset;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekg(suboffsetLoc);
  outfile.read((char *)&suboffset, sizeof(int));

  sidsize = sid.length() + 1;
  fidsize = fid.length();
  if (fidsize > 0)
    ++fidsize;

  if (numSubLoc < 0) {
    outfile.seekp(suboffset);
    outfile.write((char *)&numSub, sizeof(int));
  } else {
    outfile.seekp(numSubLoc);
    outfile.write((char *)&numSub, sizeof(int));
    outfile.seekp(suboffset);
  }

  outfile.write((char *)&sidsize, sizeof(int));
  outfile.write((char *)&fidsize, sizeof(int));
  WriteBDString(outfile, sid);
  WriteBDString(outfile, fid);

  snpoffset = outfile.tellp();
  outfile.seekp(snpoffsetLoc);
  outfile.write((char *)&snpoffset, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write a the SNP information for format 4 *********************//

// [[Rcpp::export]]
int WriteBDSNPInfoC(std::string &filename,
                    int numSNPs,
                    std::string &snpid,
                    std::string &chromosome,
                    Rcpp::IntegerVector &location,
                    std::string &reference,
                    std::string &alternate,
                    Rcpp::NumericVector &aaf,
                    Rcpp::NumericVector &maf,
                    Rcpp::NumericVector &avgCall,
                    Rcpp::NumericVector &rsq,
                    int numSNPloc,
                    int snpOptionsLoc,
                    int snpOffsetLoc,
                    int nextOffsetLoc) {
  std::fstream outfile;
  int snpoffset, nextoffset;
  int stringSize[4];
  int snpOptions = 0;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  if (snpid.length() != 0)
    snpOptions |= 0x02;
  if (chromosome.size() !=  0) {
    snpOptions |= 0x04;
    if (chromosome.find('\t') == std::string::npos)
      snpOptions |= 0x08;
  }
  if (location.length() != 0)
    snpOptions |= 0x10;
  if (reference.length() != 0)
    snpOptions |= 0x20;
  if (alternate.length() != 0)
    snpOptions |= 0x40;
  if (aaf.length() != 0)
    snpOptions |= 0x80;
  if (maf.length() != 0)
    snpOptions |= 0x100;
  if (avgCall.length() != 0)
    snpOptions |= 0x200;
  if (rsq.length() != 0)
    snpOptions |= 0x400;

  outfile.seekg(snpOffsetLoc);
  outfile.read((char *)&snpoffset, sizeof(int));

  if (snpOptionsLoc < 0) {
    outfile.seekp(snpoffset);
    outfile.write((char *)&numSNPs, sizeof(int));
    outfile.write((char *)&snpOptions, sizeof(int));
  } else {
    outfile.seekp(numSNPloc);
    outfile.write((char *)&numSNPs, sizeof(int));
    outfile.seekp(snpOptionsLoc);
    outfile.write((char *)&snpOptions, sizeof(int));
    outfile.seekp(snpoffset);
  }

  stringSize[0] = snpid.length() == 0 ? 0 : snpid.length() + 1;
  stringSize[1] = chromosome.length() == 0 ? 0 : chromosome.length() + 1;
  stringSize[2] = reference.length() == 0 ? 0 : reference.length() + 1;
  stringSize[3] = alternate.length()  == 0 ? 0 : alternate.length() + 1;
  outfile.write((char *)&stringSize[0], sizeof(stringSize));

  WriteBDString(outfile, snpid);
  WriteBDString(outfile, chromosome);
  WriteBDInteger(outfile, location);
  WriteBDString(outfile, reference);
  WriteBDString(outfile, alternate);
  WriteBDNumeric(outfile, aaf);
  WriteBDNumeric(outfile, maf);
  WriteBDNumeric(outfile, avgCall);
  WriteBDNumeric(outfile, rsq);

  nextoffset = outfile.tellp();
  outfile.seekp(nextOffsetLoc);
  outfile.write((char *)&nextoffset, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write a the index information for format 3.4 *****************//

// [[Rcpp::export]]
int WriteBDIndexArray3_4(std::string &filename,
                         int numSNPs) {
  std::fstream outfile;
  const int zero = 0;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekg(72);

  for (int i = 0; i < numSNPs; ++i)
    outfile.write((char *)&zero, sizeof(int));

  outfile.close();
  return 0;
}

//  ********** Write a the index information for format 4.4 *****************//

// [[Rcpp::export]]
int WriteBDIndexArray4(std::string &filename,
                      int numSNPs,
                      int indexoffsetLoc,
                      int dosageoffsetloc) {
  std::fstream outfile;
  int indexoffset, dosageoffset;
  const int zero = 0;

  // Open the file for appending
  // Only opens for input and output
  outfile.open(filename.c_str(), READWRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  outfile.seekg(indexoffsetLoc);
  outfile.read((char *)&indexoffset, sizeof(int));
  outfile.seekp(indexoffset);

  for (int i = 0; i < numSNPs; ++i)
    outfile.write((char *)&zero, sizeof(int));

  dosageoffset = outfile.tellp();
  outfile.seekp(dosageoffsetloc);
  outfile.write((char *)&dosageoffset, sizeof(int));

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

extern const int NUMBEROFBASES = 3;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfff3 is 65,534 or 2^32 - 1
// 0x2710 is 10,000
extern const unsigned short USBASE[NUMBEROFBASES] = {
  0x7ffe, // Used for format 1.1
  0xfffe, // Used for format 1.2
  0x2710  // Used for all other formats
};

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const double DBASE[NUMBEROFBASES] = {
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
//    if (i < 10)
//      Rcpp::Rcout << x[i] << '\t'
//                  << *ps1 << '\t'
//                  << r1 << '\t'
//                  << x1 << '\t'
//                  << r2 << '\t'
//                  << x2 << std::endl;
  }
}

// Write only the dosage data
// Paramater filename - name of bindary dosage file
// Parameter dosage - vector of dosage values to write
// Parameter usdosage - vector used to store the converted dosage values. This
//                      is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
// [[Rcpp::export]]
int WriteBinaryDosageDataC(const std::string &filename,
                          Rcpp::NumericVector &dosage,
                          Rcpp::IntegerVector &us,
                          int base) {
  std::ofstream outfile;

  // Opens file for writing binary data by appending
  outfile.open(filename.c_str(), WRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  DoubleToUShort(dosage, us, base - 1);
  outfile.write((char *)&us[0], dosage.size() * sizeof(short));

  outfile.close();
  return 0;
}

// Write only the p1 and p2 data
// Paramater filename - name of bindary dosage file
// Parameter p1 - vector of P(g=1) to write
// Parameter p2 - vector of P(g=2) to write
// Parameter us - vector used to store the converted values. This
//                is passed to avoid allocating and dealllocating memory
// Parameter base - Index of USBASE to use as base
// [[Rcpp::export]]
int WriteBinaryP1P2Data(const std::string &filename,
                            Rcpp::NumericVector &p1,
                            Rcpp::NumericVector &p2,
                            Rcpp::IntegerVector &us,
                            int base) {
  std::ofstream outfile;

  // Opens file for writing binary data by appending
  outfile.open(filename.c_str(), WRITEBINARY);
  if (!outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  DoubleToUShort(p1, us, base - 1);
  outfile.write((char *)&us[0], p1.size() * sizeof(short));
  DoubleToUShort(p2, us, base - 1);
  outfile.write((char *)&us[0], p2.size() * sizeof(short));

  outfile.close();
  return 0;
}
