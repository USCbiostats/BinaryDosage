#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

const std::ios_base::openmode READBINARY = std::ios_base::in | std::ios_base::binary;

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
  int i, j;
  int format, subformat;

  // Open the file - if file already exists, truncates to size 0.
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  infile.read((char *)&magicWord, sizeof(int));
  infile.read((char *)&formatInt, sizeof(int));
  if (magicWord != MAGICWORD)
    Rcpp::Rcout << "File does not appear to be a binary dosage file" << std::endl;

  format = 0;
  subformat = 0;
  for (i = 0; i < FORMAT.size(); ++i) {
    for (j = 0; j < FORMAT[i].size(); ++j) {
      if(FORMAT[i][j] == formatInt) {
        format = i + 1;
        subformat = j + 1;
        i = FORMAT.size();
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
int ReadBinaryDosageHeader3A(std::string &filename) {
  std::ifstream infile;
  int numSubjects = 0;

  // Open the file for reading
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  infile.seekg(8);
  infile.read((char *)&numSubjects, sizeof(int));

  infile.close();
  return numSubjects;
}

// Writes the additional header info for formats 3.3 and 3.4
// Parameter filename - Name of binary dosage file
// Return - StringVector with MD5 hashes of family and map files
// [[Rcpp::export]]
Rcpp::StringVector ReadBinaryDosageHeader3B(std::string &filename) {
  std::ifstream infile;
  Rcpp::StringVector md5(2);
  char md5hash[33];
  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return 1;
  }

  infile.seekg(8);
  infile.read(md5hash, 32);
  md5hash[32] = 0;
  md5[0] = md5hash;
  infile.read(md5hash, 32);
  md5[1] = md5hash;

  infile.close();
  return md5;
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
  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }


  infile.read((char *)&numSubjects, sizeof(int));
  infile.read((char *)&numSNPs, sizeof(int));
  infile.read((char *)&numGroups, sizeof(int));
  infile.read((char *)&subjectOptions, sizeof(int));
  infile.read((char *)&snpOptions, sizeof(int));
  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  infile.close();
  return (Rcpp::List::create(Rcpp::Named("numsub") = numSubjects,
                             Rcpp::Named("numSNPs") = numSNPs,
                             Rcpp::Named("numgroups") = numGroups,
                             Rcpp::Named("suboptions") = subjectOptions,
                             Rcpp::Named("snpOptions") = snpOptions,
                             Rcpp::Named("subjectoffset") = subjectOffset,
                             Rcpp::Named("snpoffset") = snpOffset,
                             Rcpp::Named("dosageoffset") = dosageOffset));
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

  Rcpp::List retVal;

  // Open the file for appending
  // Only opens for output
  infile.open(filename.c_str(), READBINARY);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    return retVal;
  }

  infile.read((char *)&subjectOffset, sizeof(int));
  infile.read((char *)&snpOffset, sizeof(int));
  infile.read((char *)&indexOffset, sizeof(int));
  infile.read((char *)&dosageOffset, sizeof(int));

  infile.close();
  return Rcpp::List::create(Rcpp::Named("suboffset") = subjectOffset,
                            Rcpp::Named("snpoffset") = snpOffset,
                            Rcpp::Named("indexoffset") = indexOffset,
                            Rcpp::Named("dosageoffset") = dosageOffset);
}
