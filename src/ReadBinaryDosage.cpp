#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

extern const std::ios_base::openmode READWRITEBINARY;
extern const std::ios_base::openmode READBINARY;
extern const std::ios_base::openmode WRITEBINARY;

extern const int NUMBEROFBASES;
extern const unsigned short USBASE[];
extern const double DBASE[];

//***************************************************************************//
//                        Reading the header                                 //
//***************************************************************************//
// These functions read the headers to the binary dosage files

//  ************************ Constants **************************************//
// Magic word for binary dosage files
extern const int MAGICWORD;
// Format ID stored in file
extern const std::vector<std::vector<int> > FORMAT;

// Reads the base header for a binary dosage file
// Parameter filename - Name of binary dosage file
// Return - list with format and subformat
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageBaseHeader(std::string &filename) {
  std::ifstream infile;
  int magicWord;
  int formatInt;
  unsigned int ui, uj;
  int format, subformat;
  Rcpp::List retVal;

  // Open the file - if file already exists, truncates to size 0.
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.read((char *)&magicWord, sizeof(int));
  infile.read((char *)&formatInt, sizeof(int));
  if (magicWord != MAGICWORD) {
    Rcpp::Rcerr << "File does not appear to be a binary dosage file" << std::endl;
    return retVal;
  }

  format = 0;
  subformat = 0;
  for (ui = 0; ui < FORMAT.size(); ++ui) {
    for (uj = 0; uj < FORMAT[ui].size(); ++uj) {
      if(FORMAT[ui][uj] == formatInt) {
        format = ui + 1;
        subformat = uj + 1;
        ui = FORMAT.size();
        break;
      }
    }
  }

  infile.close();
  return Rcpp::List::create(Rcpp::Named("format") = format,
                            Rcpp::Named("subformat") = subformat);
}

// Writes the additional header info for formats 3.1 and 3.2
// Parameter filename - Name of binary dosage file
// Return - number of subjects in data
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader3A(std::string &filename) {
  std::ifstream infile;
  int numSubjects = 0;
  Rcpp::List retVal;

  // Open the file for reading
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));

  infile.close();
  return Rcpp::List::create(Rcpp::Named("numsub") = numSubjects);
}

// Writes the additional header info for formats 3.3 and 3.4
// Parameter filename - Name of binary dosage file
// Return - StringVector with MD5 hashes of family and map files
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader3B(std::string &filename) {
  std::ifstream infile;
  Rcpp::StringVector md5(2);
  Rcpp::List retVal;
  char md5hash[33];

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.seekg(8);
  infile.read(md5hash, 32);
  md5hash[32] = 0;
  md5[0] = md5hash;
  infile.read(md5hash, 32);
  md5[1] = md5hash;

  infile.close();
  return Rcpp::List::create(Rcpp::Named("md5") = md5);
}

std::string ReadBDString(std::ifstream &infile, int length) {
  char *stringRead = NULL;
  std::string outString = "";

  if (length != 0) {
    stringRead = new char[length];
    infile.read(stringRead, length);
    outString = stringRead;
  }
  return outString;
}

std::vector<int> ReadBDInteger(std::ifstream &infile, int length) {
  std::vector<int> retVal;

  if (length > 0) {
    retVal.resize(length);
    infile.read((char *)retVal.data(), length * sizeof(int));
  }
  return retVal;
}

std::vector<double> ReadBDNumeric(std::ifstream &infile, int length) {
  std::vector<double> retVal;

  if (length > 0) {
    retVal.resize(length);
    infile.read((char *)retVal.data(), length * sizeof(double));
  }
  return retVal;
}

Rcpp::List ReadBDSubjects(std::ifstream &infile) {
  int sidsize;
  int fidsize;
  std::string sidString;
  std::string fidString;

  infile.read((char *)&sidsize, sizeof(int));
  infile.read((char *)&fidsize, sizeof(int));

  sidString = ReadBDString(infile, sidsize);
  fidString = ReadBDString(infile, fidsize);

  return Rcpp::List::create(Rcpp::Named("sidsize") = sidsize,
                            Rcpp::Named("fidsize") = fidsize,
                            Rcpp::Named("sidstring") = sidString,
                            Rcpp::Named("fidstring") = fidString);
}

Rcpp::List ReadBDSNPs(std::ifstream &infile, int numSNPs, int numGroups, int snpoptions) {
  int snpSize, chrSize, refSize, altSize;
  std::string snpString, chrString, refString, altString;
  std::vector<int> location;
  std::vector<double> aaf, maf, avgCall, rsq;

  infile.read((char *)&snpSize, sizeof(int));
  infile.read((char *)&chrSize, sizeof(int));
  infile.read((char *)&refSize, sizeof(int));
  infile.read((char *)&altSize, sizeof(int));

  snpString = ReadBDString(infile, snpSize);
  chrString = ReadBDString(infile, chrSize);
  if ((snpoptions & 0x0010) != 0)
    location = ReadBDInteger(infile, numSNPs);
  refString = ReadBDString(infile, refSize);
  altString = ReadBDString(infile, altSize);

  if ((snpoptions & 0x0080) != 0)
    aaf = ReadBDNumeric(infile, numSNPs * numGroups);
  if ((snpoptions & 0x0100) != 0)
    maf = ReadBDNumeric(infile, numSNPs * numGroups);
  if ((snpoptions & 0x0200) != 0)
    avgCall = ReadBDNumeric(infile, numSNPs * numGroups);
  if ((snpoptions & 0x0400) != 0)
    rsq = ReadBDNumeric(infile, numSNPs * numGroups);

  return Rcpp::List::create(Rcpp::Named("snpsize") = snpSize,
                            Rcpp::Named("chrsize") = chrSize,
                            Rcpp::Named("refsize") = refSize,
                            Rcpp::Named("altsize") = altSize,
                            Rcpp::Named("snpstring") = snpString,
                            Rcpp::Named("chrstring") = chrString,
                            Rcpp::Named("location") = location,
                            Rcpp::Named("refstring") = refString,
                            Rcpp::Named("altstring") = altString,
                            Rcpp::Named("aaf") = aaf,
                            Rcpp::Named("maf") = maf,
                            Rcpp::Named("avgcall") = avgCall,
                            Rcpp::Named("rsq") = rsq);
}

// Writes the additional header info for formats 4.1 and 4.2
// Parameter filename - Name of binary dosage file
// Return - list of values in header
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader4A(std::string &filename) {
  std::ifstream infile;
  int numSubjects;
  int numSNPs;
  int numGroups;
  int subjectOptions;
  int snpOptions;
  int subjectOffset;
  int snpOffset;
  int dosageOffset;
  std::vector<int> groupSizes;
  Rcpp::List samples;
  Rcpp::List snps;
  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));
  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&numGroups, sizeof(int));
  infile.read((char *)&subjectOptions, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  groupSizes.resize(numGroups);
  infile.read((char *)groupSizes.data(), numGroups * sizeof(int));

  samples = ReadBDSubjects(infile);
  snps = ReadBDSNPs(infile, numSNPs, numGroups, snpOptions);

  infile.close();

  return Rcpp::List::create(Rcpp::Named("numsub") = numSubjects,
                            Rcpp::Named("numSNPs") = numSNPs,
                            Rcpp::Named("numgroups") = numGroups,
                            Rcpp::Named("suboptions") = subjectOptions,
                            Rcpp::Named("snpoptions") = snpOptions,
                            Rcpp::Named("subjectoffset") = subjectOffset,
                            Rcpp::Named("snpoffset") = snpOffset,
                            Rcpp::Named("dosageoffset") = dosageOffset,
                            Rcpp::Named("groups") = groupSizes,
                            Rcpp::Named("samples") = samples,
                            Rcpp::Named("snps") = snps);
}


// Writes the additional header info for formats 4.2 and 4.3
// Parameter filename - Name of binary dosage file
// Return - 0 successful, 1 error
// [[Rcpp::export]]
Rcpp::List ReadBinaryDosageHeader4B (std::string &filename) {
  std::ifstream infile;
  int subjectOffset;
  int snpOffset;
  int indexOffset;
  int dosageOffset;
  int numGroups;
  std::vector<int> groupSizes;
  int snpOptions;
  int numSub;
  Rcpp::List samples;
  int numSNPs;
  Rcpp::List snps;
  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.seekg(8);
  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&indexOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  infile.read((char *)&numGroups, sizeof(int));
  groupSizes.resize(numGroups);
  infile.read((char *)groupSizes.data(), numGroups * sizeof(int));

  infile.read((char *)&numSub, sizeof(int));
  samples = ReadBDSubjects(infile);

  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  snps = ReadBDSNPs(infile, numSNPs, numGroups, snpOptions);

  infile.close();

  return Rcpp::List::create(Rcpp::Named("suboffset") = subjectOffset,
                            Rcpp::Named("snpoffset") = snpOffset,
                            Rcpp::Named("indexoffset") = indexOffset,
                            Rcpp::Named("dosageoffset") = dosageOffset,
                            Rcpp::Named("numgroups") = numGroups,
                            Rcpp::Named("groups") = groupSizes,
                            Rcpp::Named("numsub") = numSub,
                            Rcpp::Named("samples") = samples,
                            Rcpp::Named("numSNPs") = numSNPs,
                            Rcpp::Named("snpoptions") = snpOptions,
                            Rcpp::Named("snps") = snps);
}

// [[Rcpp::export]]
Rcpp::IntegerVector ReadBDIndicesS4(std::string filename,
                                    int numSNPs,
                                    int indexStart) {
  std::ifstream infile;
  Rcpp::IntegerVector retVal(numSNPs);

  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.seekg(indexStart);
  infile.read((char *)&retVal[0], numSNPs * sizeof(int));

  infile.close();

  return retVal;
}

// Routine to convert short values to double
// Parameter us - is the vector of shorts to convert
// Parameter x - vector of unsigned doubles to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff for shorts
void UShortToDouble(Rcpp::IntegerVector &us,
                    Rcpp::NumericVector &x,
                    const int base) {
  unsigned short *ps1;
  double *d;
  int i;

  ps1 = (unsigned short *)&us[0];
  d = (double *)&x[0];
  for (i = 0; i < x.size(); ++i, ++ps1, ++d) {
    if (*ps1 == 0xffff) {
      // Missing
      *d = NA_REAL;
    } else {
      *d = DBASE[base] * *ps1;
    }
  }
}

// [[Rcpp::export]]
int ReadBinaryDosageDataC(std::string &filename,
                          int headersize,
                          int snp,
                          Rcpp::NumericVector &dosage,
                          Rcpp::IntegerVector &usdosage,
                          int base) {
  std::ifstream infile;
  std::streampos loc;

  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open input file" << std::endl;
    return 1;
  }

  loc = headersize + 2 * (snp - 1) * dosage.size();
  infile.seekg(loc);
  infile.read((char *)&usdosage[0], usdosage.size() * sizeof(short));
  UShortToDouble(usdosage, dosage, base - 1);
  infile.close();
  return 0;
}

// [[Rcpp::export]]
int ReadBinaryDosageDataP1P2(std::string &filename,
                             int headersize,
                             int snp,
                             Rcpp::NumericVector &dosage,
                             Rcpp::NumericVector &p0,
                             Rcpp::NumericVector &p1,
                             Rcpp::NumericVector &p2,
                             Rcpp::IntegerVector &us,
                             int base) {
  std::ifstream infile;
  std::streampos loc;

  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open input file" << std::endl;
    return 1;
  }

  loc = headersize + 4 * (snp - 1) * dosage.size();
  infile.seekg(loc);
  infile.read((char *)&us[0], us.size() * sizeof(short));
  UShortToDouble(us, p1, base - 1);
  infile.read((char *)&us[0], us.size() * sizeof(short));
  UShortToDouble(us, p2, base - 1);
  dosage = p1 + p2 + p2;
  p0 = 1. - p1 - p2;
  infile.close();
  return 0;
}
