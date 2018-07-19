#include <RcppArmadillo.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "BinaryDosageUtilities.h"
#include "BinaryDosageWrite.h"
#include "VCFtoBinaryDosage.h"

const char CReadBinaryDosage::m_magicWord42[8] = { 'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x2};

// ********************************************************************************************
//
// Reading in the binary dosage file information
//
// ********************************************************************************************
CReadBinaryDosage::CReadBinaryDosage() {
  m_numSubjects = 0;
  m_numSNPs = 0;
  m_numGroups = 0;
  m_subjectOptions = 0;
  m_SNPOptions = 0;
  m_subjectStart = 0;
  m_SNPStart = 0;
  m_SNPNameSize = 0;
  m_chromosomeSize = 0;
  m_refAlleleSize = 0;
  m_altAlleleSize = 0;
  m_dosageStart = 0;
  m_subjectIDSize = 0;
  m_familyIDSize = 0;
}

CReadBinaryDosage::~CReadBinaryDosage() {
  m_infile.close();
}

int CReadBinaryDosage::ReadHeader() {
  m_infile.read(m_magicWord, 8);
  if (std::memcmp(m_magicWord, m_magicWord42, 8)) {
    Rcpp::Rcerr << "File does not appear to be a binary dosage file" << std::endl;
    return 1;
  }
  m_infile.read((char *)&m_numSubjects, sizeof(unsigned int));
  m_infile.read((char *)&m_numSNPs, sizeof(unsigned int));
  m_infile.read((char *)&m_numGroups, sizeof(unsigned int));
  m_infile.read((char *)&m_subjectOptions, sizeof(unsigned int));
  m_infile.read((char *)&m_SNPOptions, sizeof(unsigned int));
  m_infile.read((char *)&m_subjectStart, sizeof(unsigned int));
  m_infile.read((char *)&m_SNPStart, sizeof(unsigned int));
  m_infile.read((char *)&m_dosageStart, sizeof(unsigned int));
  std::vector<unsigned int>().swap(m_groupSize);
  m_groupSize.resize(m_numGroups);
  m_infile.read((char *)&m_groupSize[0], m_numGroups * sizeof(unsigned int));

//  Rcpp::Rcout << m_numSubjects << '\t' << m_numSNPs << '\t' << m_numGroups << std::endl;
//  Rcpp::Rcout << std::hex << m_subjectOptions << '\t' << m_SNPOptions << '\n' << std::dec << m_subjectStart
//              << '\t' << m_SNPStart << '\t' << m_dosageStart << std::endl;

  if (m_infile.fail())
    return 1;

  return 0;
}

int CReadBinaryDosage::ReadSubjects() {
  std::vector<char> idArray;
  std::istringstream iss;
  unsigned int ui;

  std::vector<std::string>().swap(m_subjectID);
  std::vector<std::string>().swap(m_familyID);

  m_infile.seekg(m_subjectStart);
  m_infile.read((char *)&m_subjectIDSize, sizeof(unsigned int));
  m_infile.read((char *)&m_familyIDSize, sizeof(unsigned int));

  idArray.resize(m_subjectIDSize);
  m_infile.read(&idArray[0], m_subjectIDSize);
  m_subjectID.resize(m_numSubjects);
  iss.str(&idArray[0]);
  for (ui = 0; ui < m_numSubjects; ++ui)
    iss >> m_subjectID[ui];

  if (m_familyIDSize > 0) {
    std::vector<char>().swap(idArray);
    idArray.resize(m_familyIDSize);
    iss.str(&idArray[0]);
    iss.clear();
    for (ui = 0; ui < m_numSubjects; ++ui)
      iss >> m_familyID[ui];
  }
//  Rcpp::Rcout << m_infile.tellg() << std::endl;

  return 0;
}

int CReadBinaryDosage::ReadSNPInfo() {
  std::vector<char> idArray;
  std::istringstream iss;
  std::string rchr;
  unsigned int ui;

  std::vector<std::string>().swap(m_SNPName);
  std::vector<std::string>().swap(m_chromosome);
  std::vector<std::string>().swap(m_refAllele);
  std::vector<std::string>().swap(m_altAllele);
  std::vector<unsigned int>().swap(m_location);
  std::vector<double>().swap(m_altFreq);
  std::vector<double>().swap(m_maf);
  std::vector<double>().swap(m_avgCall);
  std::vector<double>().swap(m_rSquared);

  m_infile.seekg(m_SNPStart);
  m_infile.read((char *)&m_SNPNameSize, sizeof(unsigned int));
  m_infile.read((char *)&m_chromosomeSize, sizeof(unsigned int));
  m_infile.read((char *)&m_refAlleleSize, sizeof(unsigned int));
  m_infile.read((char *)&m_altAlleleSize, sizeof(unsigned int));

  if (m_SNPNameSize > 0) {
    std::vector<char>().swap(idArray);
    idArray.resize(m_SNPNameSize);
    m_infile.read(&idArray[0], m_SNPNameSize);
    iss.str(&idArray[0]);
    iss.clear();
    m_SNPName.resize(m_numSNPs);
    for (ui = 0; ui < m_numSNPs; ++ui)
      iss >> m_SNPName[ui];
  }

  if (m_chromosomeSize > 0) {
    std::vector<char>().swap(idArray);
    idArray.resize(m_chromosomeSize);
    m_infile.read(&idArray[0], m_chromosomeSize);
    iss.str(&idArray[0]);
    iss.clear();
    m_chromosome.resize(m_numSNPs);
    if (m_SNPOptions & 0x0008) {
      iss >> rchr;
      for (ui = 0; ui < m_numSNPs; ++ui)
        m_chromosome[ui] = rchr;
    } else {
      for (ui = 0; ui < m_numSNPs; ++ui)
        iss >> m_chromosome[ui];
    }
  }

  m_location.resize(m_numSNPs);
  m_infile.read((char *)&m_location[0], m_numSNPs * sizeof(unsigned int));

  if (m_refAlleleSize > 0) {
    std::vector<char>().swap(idArray);
    idArray.resize(m_refAlleleSize);
    m_infile.read(&idArray[0], m_refAlleleSize);
    iss.str(&idArray[0]);
    iss.clear();
    m_refAllele.resize(m_numSNPs);
    for (ui = 0; ui < m_numSNPs; ++ui)
      iss >> m_refAllele[ui];
  }

  if (m_altAlleleSize > 0) {
    std::vector<char>().swap(idArray);
    idArray.resize(m_altAlleleSize);
    m_infile.read(&idArray[0], m_altAlleleSize);
    iss.str(&idArray[0]);
    iss.clear();
    m_altAllele.resize(m_numSNPs);
    for (ui = 0; ui < m_numSNPs; ++ui)
      iss >> m_altAllele[ui];
  }

  m_altFreq.resize(m_numSNPs);
  m_infile.read((char *)&m_altFreq[0], m_numSNPs * sizeof(double));
  m_maf.resize(m_numSNPs);
  m_infile.read((char *)&m_maf[0], m_numSNPs * sizeof(double));
  m_avgCall.resize(m_numSNPs);
  m_infile.read((char *)&m_avgCall[0], m_numSNPs * sizeof(double));
  m_rSquared.resize(m_numSNPs);
  m_infile.read((char *)&m_rSquared[0], m_numSNPs * sizeof(double));

//  Rcpp::Rcout << m_infile.tellg() << std::endl;

  return 0;
}

// Read the the amount of memory used for each SNP dosase information
int CReadBinaryDosage::ReadSNPDataSize() {
  unsigned int ui;

  std::vector<unsigned int>().swap(m_SNPDataSize);
  m_SNPDataSize.resize(m_numSNPs);

  m_infile.seekg(m_dosageStart);
  for (ui = 0; ui < m_numSNPs; ++ui) {
    m_infile.read((char *)&m_SNPDataSize[ui], sizeof(unsigned int));
    m_infile.seekg(m_SNPDataSize[ui], std::ios::cur);
  }
  return 0;
}

int CReadBinaryDosage::ReadFileInfo(const std::string &binaryDosageFilename) {
  m_infile.close();
  m_infile.clear();

  m_infile.open(binaryDosageFilename.c_str(), std::ios::in | std::ios::binary);
  if (!m_infile.good()) {
    Rcpp::Rcerr << "Unable to open file:\t" << binaryDosageFilename << std::endl;
    return 1;
  }

  if (ReadHeader())
    return 1;
  if (ReadSubjects())
    return 1;
  if (ReadSNPInfo())
    return 1;
  if (ReadSNPDataSize())
    return 1;

  return 0;
}


//' Function to convert a VCF file to a binary dosage file for GxEScan
//'
//' Function to convert a VCF file to a binary dosage file for GxEScan.
//' The binary formatted file is smaller that the VCF file and reads in much quicker.
//'
//' @param vcfFilename
//' Name of VCF file
//' @param infoFilename
//' Name of information file
//' @param outputFilename
//' Name of the binary dosage file
//' @return
//' 0 Success
//' 1 Failure
//' @export
// [[Rcpp::export]]
int VCF2BD_C(const std::string &vcfFilename, const std::string &infoFilename, const std::string &outputFilename) {
  CWriteBinaryDosageFromVCF wvcf;

  return wvcf.WriteBinaryDosageFile(vcfFilename, infoFilename, outputFilename);
}

//' Function to convert a VCF file from HRC to a binary dosage file for GxEScan without an info file
//'
//' Function to convert a VCF file from HRC to a binary dosage file for GxEScan without an info file.
//' The binary formatted file is smaller that the VCF file and reads in much quicker.
//'
//' @param vcfFilename
//' Name of VCF file
//' @param outputFilename
//' Name of the binary dosage file
//' @return
//' 0 Success
//' 1 Failure
//' @export
// [[Rcpp::export]]
int VCF53toBD_C(const std::string &vcfFilename, const std::string &bdFilename) {
  CVCFtoBinaryDosage42 vcfbd42;

  vcfbd42.Convert(vcfFilename, bdFilename);
  return 0;
}
//' Function to convert an Impute 2 dosage file to a binary dosage file for GxEScan
//'
//' Function to convert a VCF file to a binary dosage file for GxEScan.
//' The Impute 2 dosage file has to be created from a VCF file using VCF tools
//' The binary formatted file is smaller that the VCF file and reads in much quicker.
//'
//' @param vcfFilename
//' Name of VCF file
//' @param infoFilename
//' Name of information file
//' @param smapleFilename
//' Name of sample file - file with subject IDs
//' @param outputFilename
//' Name of the binary dosage file
//' @return
//' 0 Success
//' 1 Failure
//' @export
// [[Rcpp::export]]
int VCF2Impute2BD(const std::string &vcfFilename, const std::string &infoFilename, const std::string sampleFilename, const std::string &outputFilename) {
  CWriteBinaryDosageFromVCF2Impute wvcf;

  return wvcf.WriteBinaryDosageFile(vcfFilename, infoFilename, sampleFilename, outputFilename);
}

//' Function to convert a VCF file to a binary dosage file for GxEScan
//'
//' Function to convert a VCF file to a binary dosage file for GxEScan.
//' The binary formatted file is smaller that the VCF file and reads in much quicker.
//'
//' @param binaryDosagFilename
//' Name of binary dosage file
//' @return
//' 0 Success
//' 1 Failure
//' @export
// [[Rcpp::export]]
Rcpp::List ReadBDInfo(const std::string &binaryDosageFilename) {
  Rcpp::List bdData;
  Rcpp::List subjectList;
  Rcpp::List snpList;
  Rcpp::List fileData;
  CReadBinaryDosage rbd;
  std::vector<std::string> names;

  if (rbd.ReadFileInfo(binaryDosageFilename))
    return bdData;

  if (rbd.FamilyID().size() == 0) {
    subjectList = Rcpp::List::create(Rcpp::Named("IID") = rbd.SubjectID());
    names.resize(1);
    names[0] = "IID";
  } else {
    subjectList = Rcpp::List::create(Rcpp::Named("FID") = rbd.FamilyID(), Rcpp::Named("IID") = rbd.SubjectID());
    names.resize(2);
    names[0] = "FID";
    names[1] = "IID";
  }

  Rcpp::DataFrame subjectDF(subjectList);
  subjectDF.attr("names") = names;

  fileData = Rcpp::List::create(Rcpp::Named("FileName") = binaryDosageFilename,
                                Rcpp::Named("NumSubject") = rbd.NumSubjects(),
                                Rcpp::Named("NumSNPs") = rbd.NumSNPs(),
                                Rcpp::Named("NumGroups") = rbd.NumGroups(),
                                Rcpp::Named("SubjectOptions") = rbd.SubjectOptions(),
                                Rcpp::Named("SNPOptions") = rbd.SNPOptions(),
                                Rcpp::Named("GroupSize") = rbd.GroupSize(),
                                Rcpp::Named("SubjectStart") = rbd.SubjectStart(),
                                Rcpp::Named("SNPStart") = rbd.SNPStart(),
                                Rcpp::Named("DosageStart") = rbd.DosageStart(),
                                Rcpp::Named("DataSize") = rbd.SNPDataSize());
  snpList = Rcpp::List::create(Rcpp::Named("Chromosome") = rbd.Chromosome(),
                               Rcpp::Named("Location") = rbd.Location(),
                               Rcpp::Named("RefAllele") = rbd.RefAllele(),
                               Rcpp::Named("AltAllele") = rbd.AltAllele(),
                               Rcpp::Named("AltFreq") = rbd.AltFreq(),
                               Rcpp::Named("maf") = rbd.MAF(),
                               Rcpp::Named("AvgCall") = rbd.AvgCall(),
                               Rcpp::Named("rSquared") = rbd.RSquared());

  Rcpp::DataFrame SNPDF(snpList);

  bdData = Rcpp::List::create(Rcpp::Named("Subjects") = subjectDF,
                              Rcpp::Named("SNPs") = SNPDF,
                              Rcpp::Named("FileData") = fileData);

  return bdData;
}

//' Function read binary dosage file SNPs into a data.table
//'
//' Function read binary dosage file SNPs into a data.table
//'
//' @param filename
//' Name of binary dosage file
//' @param numSubjects
//' Number of subjects in file
//' @param numSNPs
//' Number of SNPs in file
//' @param dosageStart
//' Where dosage values start in file
//' @param dataSize
//' Size of SNP data on disk
//' @param snps
//' List of SNPs to extract
//' @param dosageptr
//' Pointers to dosage values in data.table
//' @param p0ptr
//' Pointers to p0 values in data.table
//' @param p1ptr
//' Pointers to p1 values in data.table
//' @param p2ptr
//' Pointers to p2 values in data.table
//' @return
//' 0 Success
//' 1 failure
//' @export
// [[Rcpp::export]]
int ReadBDSNPs(std::string &filename, unsigned int numSubjects, unsigned int numSNPs,
               unsigned int dosageStart, Rcpp::IntegerVector &dataSize, Rcpp::IntegerVector &snps,
               Rcpp::NumericVector &dosageptr, Rcpp::NumericVector &p0ptr, Rcpp::NumericVector &p1ptr, Rcpp::NumericVector &p2ptr) {
  std::vector<short> ud(4*numSubjects);
  unsigned int ev;
  std::ifstream infile;
  unsigned int ui, uj, uk;
  unsigned int arraySize;
  double *d;
  double *p0, *p1, *p2;
  std::ios::streampos snpPos;

  infile.open(filename.c_str(), std::ios::in | std::ios::binary);
  if (!infile.good()) {
    Rcpp::Rcerr << "Unable to open file:\t" << filename << std::endl;
    return 1;
  }

  for (ui = 0; ui < snps.size(); ++ui) {
    if (ui == 0 || snps[ui - 1] >= snps[ui]) {
      if (ui != 0)
        Rcpp::Rcout << "Resetting\t" << snps[ui - 1] << '\t' << snps[ui] << std::endl;
      infile.seekg(dosageStart);
      uk = 1;
    } else {
      ++uk;
    }

    memcpy(&d, &dosageptr[ui], sizeof(double *));
    memcpy(&p0, &p0ptr[ui], sizeof(double *));
    memcpy(&p1, &p1ptr[ui], sizeof(double *));
    memcpy(&p2, &p2ptr[ui], sizeof(double *));
//    Rcpp::Rcout << d << '\t' << p0 << '\t' << p1 << '\t' << p2 << std::endl;
//    snpPos = infile.tellg();
    snpPos = 0;
//    Rcpp::Rcout << infile.tellg() << '\t';
    for (; uk < (unsigned int)snps[ui]; ++uk) {
//      infile.read((char *)&arraySize, sizeof(unsigned int));
//      infile.seekg(arraySize, std::ios::cur);
      snpPos += dataSize[uk - 1] + 4;
//      Rcpp::Rcout << "Skipped array size\t" << arraySize << std::endl;
    }
//    Rcpp::Rcout << snpPos << '\t' << infile.tellg();
    infile.seekg(snpPos, std::ios::cur);
//    Rcpp::Rcout << '\t' << snpPos << '\t' << infile.tellg() << std::endl;
    infile.read((char *)&arraySize, sizeof(unsigned int));
//    Rcpp::Rcout << "Array size\t" << arraySize << std::endl;

    if (arraySize == 0) {
      for (uj = 0; uj < numSubjects; ++uj)
        d[uj] = NA_REAL;
      continue;
    }
    infile.read((char *)&ud[0], arraySize);
    ev = numSubjects;
    for (uj = 0; uj < numSubjects; ++uj) {
      if (ud[uj] == 20001) {
        d[uj] = NA_REAL;
        continue;
      }
      if (ud[uj] & 0x8000) {
        ud[uj] &= 0x7fff;
        d[uj] = ud[uj] / 10000.;
        if (ud[ev] & 0x8000) {
          ud[ev] &= 0x7fff;
          p1[uj] = ud[ev] / 10000.;
          ++ev;
          p0[uj] = ud[ev] / 10000.;
          ++ev;
          p2[uj] = ud[ev] / 10000.;
          ++ev;
        } else {
          p1[uj] = ud[ev] / 10000.;
          ++ev;
          p2[uj] = (d[uj] - p1[uj]) / 2;
          p0[uj] = 1. - p2[uj] - p1[uj];
        }
      } else {
        d[uj] = ud[uj] / 10000.;
        if (d[uj] > 1) {
          p0[uj] = 0;
          p2[uj] = d[uj] - 1;
          p1[uj] = 1 - p2[uj];
        } else {
          p2[uj] = 0;
          p1[uj] = d[uj];
          p0[uj] = 1 - p1[uj];
        }
      }
    }
//    Rcpp::Rcout << "Read SNP" << std::endl;
  }

  infile.close();
  return 0;
}

//' Function to find the pointer to a column in a data.table
//'
//' Function to find the pointer to a column in a data.table
//'
//' @param x
//' Column in data.table
//' @param y
//' Array for storing pointers
//' @param n
//' Location in array to store pointer
//' @return
//' 0 Success
//' 1 failure
//' @export
// [[Rcpp::export]]
int FindPointer(Rcpp::NumericVector &x, Rcpp::NumericVector &y, unsigned int n) {
  double *dp;

  dp = &x[0];
//  Rcpp::Rcout << dp << std::endl;
//  x[0] = 0;
  memcpy(&y[n - 1], &dp, sizeof(double *));
  return 0;
}

//' Function to print pointers from an array
//'
//' Function to print pointers from an array
//'
//' @param y
//' Array of stored pointers
//' @return
//' 0 Success
//' 1 failure
//' @export
// [[Rcpp::export]]
int PrintPointer(Rcpp::NumericVector &y) {
  double *dp;
  unsigned int ui;

  for (ui = 0; ui < y.size(); ++ui) {
    memcpy(&dp, &y[ui], sizeof(double));
    Rcpp::Rcout << dp << std::endl;
  }
  return 0;
}
