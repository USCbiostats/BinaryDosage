#include <RcppArmadillo.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "BinaryDosageUtilities.h"

const char CWriteBinaryDosageFromVCF::m_magicWord[8] = { 'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x2};
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

// ********************************************************************************************
//
// Converting VCF to binary dosage file
//
// ********************************************************************************************
CWriteBinaryDosageFromVCF::CWriteBinaryDosageFromVCF() {
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

CWriteBinaryDosageFromVCF::~CWriteBinaryDosageFromVCF() {
  // Close the files - not much else to do
  m_infile.close();
  m_outfile.close();
}

// Read in the information file - Little to no error checking
int CWriteBinaryDosageFromVCF::ReadInfoFile(const std::string &infoFilename) {
  std::string header;
  std::string snpN;
  std::string ra;
  std::string aa;
  std::string junk;
  double af, ma, ac, rs;

  // Just in case something is open
  m_infile.close();
  m_infile.clear();
  // Set SNP options for VCF format
  m_SNPOptions = 0x07f4;
  // Clear out strings
  std::vector<std::string>().swap(m_SNPName);
  std::vector<std::string>().swap(m_chromosome);
  std::vector<unsigned int>().swap(m_location);
  std::vector<std::string>().swap(m_altAllele);
  std::vector<std::string>().swap(m_refAllele);
  std::vector<double>().swap(m_altFreq);
  std::vector<double>().swap(m_maf);
  std::vector<double>().swap(m_avgCall);
  std::vector<double>().swap(m_rSquared);
  m_SNPNameSize = 0;
  m_chromosomeSize = 0;
  m_altAlleleSize = 0;
  m_refAlleleSize = 0;
  // Open the information file
  m_infile.open(infoFilename.c_str());
  // Did the file open?
  if (!m_infile.good()) {
    Rcpp::Rcerr << "Unable to open information file:\t" << infoFilename << std::endl;
    return 1;
  }
  // Read in the header line
  getline(m_infile, header);
  // Read in the first line
  m_infile >> snpN >> ra >> aa >> af >> ma;
  m_infile >> std::ws;
  // Check if average call(ac), and r-squared(rs) are missing
  if (m_infile.peek() == '-') {
    ac = NA_REAL;
    rs = NA_REAL;
  } else {
    m_infile >> ac >> rs;
  }
  // Read the rest of the line - data not used
  getline(m_infile, junk);
  while (!m_infile.fail()) {
    m_SNPName.push_back(snpN);
    m_SNPNameSize += snpN.length();
    m_refAllele.push_back(ra);
    m_refAlleleSize += ra.length();
    m_altAllele.push_back(aa);
    m_altAlleleSize += aa.length();
    m_altFreq.push_back(af);
    m_maf.push_back(ma);
    m_avgCall.push_back(ac);
    m_rSquared.push_back(rs);
    m_infile >> snpN >> ra >> aa >> af >> ma;
    m_infile >> std::ws;
    if (m_infile.peek() == '-') {
      ac = NA_REAL;
      rs = NA_REAL;
    } else {
      m_infile >> ac >> rs;
    }
    getline(m_infile, junk);
  }
  // Save some usefull information about memory usage
  m_numSNPs = m_SNPName.size();
  m_SNPNameSize += m_numSNPs;
  m_altAlleleSize += m_numSNPs;
  m_refAlleleSize += m_numSNPs;

//  Rcpp::Rcout << "Number of SNPs\t" << m_numSNPs << std::endl;
  m_infile.close();
  // Convert SNP names to chromosome and location - This is the HRC VCF format
  return GetChromosomesAndLocations();
}

int CWriteBinaryDosageFromVCF::GetChromosomesAndLocations() {
  unsigned int ui;
  std::istringstream iss;
  std::string rchr;
  unsigned rloc;
  bool oneChromosome = true;

  m_chromosome.reserve(m_numSNPs);
  m_location.reserve(m_numSNPs);

  for (ui = 0; ui < m_numSNPs; ++ui) {
    iss.clear();
    std::replace(m_SNPName[ui].begin(), m_SNPName[ui].end(), ':', '\t');
    iss.str(m_SNPName[ui]);
    iss >> rchr >> rloc;
    if (iss.fail()) {
      Rcpp::Rcerr << "SNP name is not in <chr>:<location> format" << std::endl;
      return 1;
    }
    m_chromosome.push_back(rchr);
    m_chromosomeSize += rchr.length();
    m_location.push_back(rloc);
    if (rchr != m_chromosome[0])
      oneChromosome = false;
  }
  // Check for one chromosome - set option flag and reduce array to one value
  if (oneChromosome) {
    m_SNPOptions |= 0x0008;
    std::vector<std::string>().swap(m_chromosome);
    m_chromosome.push_back(rchr);
    m_chromosomeSize = rchr.length() + 1;
  } else {
    m_chromosomeSize += m_numSNPs;
  }
  // Clear out SNPName vector - no longer needed
  std::vector<std::string>().swap(m_SNPName);
  m_SNPNameSize = 0;
  return 0;
}

int CWriteBinaryDosageFromVCF::ReadSubjectInfo() {
  std::string readLine;
  std::string columnName;
  unsigned int ui;
  std::istringstream iss;
  std::string riid;

  // Set subject options for VCF format - no family IDs
  m_subjectOptions = 0;
  // Clear out the subject and family IDs
  m_subjectIDSize = 0;
  m_familyIDSize = 0;
  std::vector<std::string>().swap(m_subjectID);
  std::vector<std::string>().swap(m_familyID);

  // The twelfth line should have the subject IDs
  for (ui = 0; ui < 12; ++ui)
    getline(m_infile, readLine);
  if (readLine[0] != '#' || readLine[1] != 'C') {
    Rcpp::Rcerr << "Twelfth line of the vcf file does not have the subject IDs" << std::endl;
    return 1;
  }

  // Skip the leading # mark
  iss.str(readLine.substr(1));
  // Read in the column names - currently not used
  for (ui = 0; ui < 9; ++ui)
    iss >> columnName;
  if (iss.fail())
    return 1;

  iss >> riid;
  while (!iss.fail()) {
    m_subjectID.push_back(riid);
    m_subjectIDSize += riid.length();
    iss >> riid;
  }

  m_numSubjects = m_subjectID.size();
  m_subjectIDSize += m_numSubjects;
//  Rcpp::Rcout << "Number of subjects:\t" << m_numSubjects << std::endl;

  return 0;
}

// Write header - stats about the file
void CWriteBinaryDosageFromVCF::WriteHeader() {
  m_outfile.write(m_magicWord, 8);
  m_outfile.write((char *)&m_numSubjects, sizeof(unsigned int));
  m_outfile.write((char *)&m_numSNPs, sizeof(unsigned int));
  m_outfile.write((char *)&m_numGroups, sizeof(unsigned int));
  m_outfile.write((char *)&m_subjectOptions, sizeof(unsigned int));
  m_outfile.write((char *)&m_SNPOptions, sizeof(unsigned int));
  m_subjectStart = 40 + m_numGroups * sizeof(unsigned int);
//  Rcpp::Rcout << "Subject start:\t" << m_subjectStart << std::endl;
  m_SNPStart = m_subjectStart + m_subjectIDSize + 2*sizeof(unsigned int);
//  Rcpp::Rcout << "SNP start:\t" << m_SNPStart << std::endl;
  m_dosageStart = m_SNPStart + 4*sizeof(unsigned int) + m_SNPNameSize + m_chromosomeSize + m_refAlleleSize + m_altAlleleSize
    + m_numSNPs * sizeof(unsigned int) + 4 * m_numSNPs * sizeof(double);
//  Rcpp::Rcout << "Dosage start:\t" << m_dosageStart << std::endl;
  m_outfile.write((char *)&m_subjectStart, sizeof(unsigned int));
  m_outfile.write((char *)&m_SNPStart, sizeof(unsigned int));
  m_outfile.write((char *)&m_dosageStart, sizeof(unsigned int));
  m_outfile.write((char *)&m_numSubjects, sizeof(unsigned int));
}

// Read and write the dosage values
int CWriteBinaryDosageFromVCF::WriteDosages() {
  const short addVal[3][10] = {{0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000}, {0, 100, 200, 300, 400, 500, 600, 700, 800, 900}, { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90}};
  bool oneChr;
  unsigned int numSNPs;
  unsigned int numMissing;
  bool atEOF;
  std::string rChr, rSNPName, rRefAllele, rAltAllele, junk, readValues;
  bool missingSNP;
  unsigned int rLoc;
  unsigned int ui, uj, uk, um;
  unsigned int additionalInfoSize;
  unsigned int dosageSize;
  unsigned int extra1, extra2;
  std::vector<short> sdosage;
  std::vector<short> extraValue;
  short p[4];
  const unsigned int uiZero = 0;
  int pSum, pDose;
  char *aPtr;

  sdosage.resize(m_numSubjects);
  extraValue.resize(3*m_numSubjects);

  oneChr = (m_chromosome.size() == 1);
  numSNPs = 0;
  numMissing = 0;
  atEOF = false;
  m_infile >> rChr >> rLoc >> rSNPName >> rRefAllele >> rAltAllele >> junk >> junk >> junk >> junk;
  m_infile >> std::ws;
  getline(m_infile, readValues);
  for (ui = 0; ui < m_numSNPs; ++ui) {
    if (!atEOF) {
      if (oneChr)
        missingSNP = !(rChr == m_chromosome[0] && rLoc == m_location[ui] && rRefAllele == m_refAllele[ui] && rAltAllele == m_altAllele[ui]);
      else
        missingSNP = !(rChr == m_chromosome[ui] && rLoc == m_location[ui] && rRefAllele == m_refAllele[ui] && rAltAllele == m_altAllele[ui]);
      if (missingSNP) {
        m_outfile.write((char *)&uiZero, sizeof(unsigned int));
        ++numMissing;
      } else {
        aPtr = &readValues[4];
        additionalInfoSize = 0;
        extra1 = 0;
        extra2 = 0;
        for (uj = 0; uj < m_numSubjects; ++uj, aPtr += 4) {
          memset(p, 0, sizeof(p));
          for (uk = 0; uk < 4; ++uk, ++aPtr) {
            switch(*aPtr) {
            case '0':
              break;
            case '1':
              p[uk] = 10000;
              break;
            case'2':
              p[uk] = 20000;
              break;
            case '.':
              if (uk == 0) {
                p[uk] = 20001;
              } else {
                Rcpp::Rcerr << "VCF format error" << std::endl;
                return 1;
              }
              break;
            default:
              Rcpp::Rcerr << "VCF format error" << std::endl;
            return 1;
            }
            ++aPtr;
            ++aPtr;
            if (p[uk] == 20001) {
              ++aPtr;
              break;
            }
            for (um = 0; um < 3; ++um, ++aPtr) {
              switch (*aPtr) {
              case '0':
                break;
              case '1':
              case '2':
              case '3':
              case '4':
              case '5':
              case '6':
              case '7':
              case '8':
              case '9':
                p[uk] += addVal[um][*aPtr - '0'];
                break;
              default:
                Rcpp::Rcerr << "VCF format error" << std::endl;
              return 1;
              }
            }
          }
          if (p[0] == 20001)
            continue;
          pSum = p[1] + p[2] + p[3];
          pDose = p[2] + p[3] + p[3];
          sdosage[uj] = p[0];
          if (pDose != p[0] || pSum != 10000) {
            sdosage[uj] |= 0x8000;
            extraValue[additionalInfoSize] = p[2] | 0x8000;
            ++additionalInfoSize;
            extraValue[additionalInfoSize] = p[1];
            ++additionalInfoSize;
            extraValue[additionalInfoSize] = p[3];
            ++additionalInfoSize;
            ++extra2;
          } else if (p[1] != 0 && p[3] != 0) {
            //            if (ui == 2)
            //              Rcpp::Rcout << p[0] << '\t' << p[1] << '\t' << p[2] << '\t' << p[3] << '\t' << pDose << '\t' << pSum << std::endl;
            sdosage[uj] |= 0x8000;
            extraValue[additionalInfoSize] = p[2];
            ++additionalInfoSize;
            ++extra1;
          }
        }
//        Rcpp::Rcout << extra1 << '\t' << extra2 << std::endl;
        ++numSNPs;
        dosageSize = m_numSubjects + m_numSubjects + additionalInfoSize + additionalInfoSize;
        m_outfile.write((char *)&dosageSize, sizeof(unsigned int));
        m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
        if (additionalInfoSize > 0)
          m_outfile.write((char *)&extraValue[0], additionalInfoSize + additionalInfoSize);
        m_infile >> rChr;
        if (m_infile.fail()) {
          if (m_infile.eof()) {
            atEOF = true;
          } else {
            Rcpp::Rcerr << "Unknown error when reading VCF file" << std::endl;
            return 1;
          }
        } else {
          m_infile >> rLoc >> rSNPName >> rRefAllele >> rAltAllele >> junk >> junk >> junk >> junk;
          m_infile >> std::ws;
          if (m_infile.fail()) {
            Rcpp::Rcerr << "Error reading SNP information in VCF file" << std::endl;
            return 1;
          } else {
            getline(m_infile, readValues);
          }
        }
      }
    } else {
      m_outfile.write((char *)&uiZero, sizeof(unsigned int));
      ++numMissing;
    }
  }
//  Rcpp::Rcout << "Number of SNPs\t" << numSNPs << std::endl << "Number missing from VCF file\t" << numMissing << std::endl;
  if (!atEOF) {
    Rcpp::Rcerr << "Did not reached EOF" << std::endl;
    Rcpp::Rcerr << "Possible extra SNPs or blank lines at end of VCF file" << std::endl;
    return 1;
  }
  return 0;
}
// Write the subject IDs (no family IDs in VCF files)
void CWriteBinaryDosageFromVCF::WriteSubjects() {
  unsigned int ui;
  const char tabChar = '\t';

  m_outfile.write((char *)&m_subjectIDSize, sizeof(unsigned int));
  m_outfile.write((char *)&m_familyIDSize, sizeof(unsigned int));
  for (ui = 0; ui < m_numSubjects; ++ui) {
    m_outfile.write(m_subjectID[ui].c_str(), m_subjectID[ui].length());
    m_outfile.write(&tabChar, 1);
  }
}

// Write SNP information - SNP Name, chromosome, location, refAllele, altAllele, altFreq, maf, avgCall, rsq
void CWriteBinaryDosageFromVCF::WriteSNPInfo() {
  unsigned int ui;
  const char tabChar = '\t';

  m_outfile.write((char *)&m_SNPNameSize, sizeof(unsigned int));
  m_outfile.write((char *)&m_chromosomeSize, sizeof(unsigned int));
  m_outfile.write((char *)&m_refAlleleSize, sizeof(unsigned int));
  m_outfile.write((char *)&m_altAlleleSize, sizeof(unsigned int));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  for (ui = 0; ui < m_SNPName.size(); ++ui) {
    m_outfile.write(m_SNPName[ui].c_str(), m_SNPName[ui].length());
    m_outfile.write(&tabChar, 1);
  }
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  for (ui = 0; ui < m_chromosome.size(); ++ui) {
    m_outfile.write(m_chromosome[ui].c_str(), m_chromosome[ui].length());
    m_outfile.write(&tabChar, 1);
  }
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  m_outfile.write((char *)&m_location[0], m_location.size() * sizeof(unsigned int));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  for (ui = 0; ui < m_refAllele.size(); ++ui) {
    m_outfile.write(m_refAllele[ui].c_str(), m_refAllele[ui].length());
    m_outfile.write(&tabChar, 1);
  }
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  for (ui = 0; ui < m_altAllele.size(); ++ui) {
    m_outfile.write(m_altAllele[ui].c_str(), m_altAllele[ui].length());
    m_outfile.write(&tabChar, 1);
  }
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  m_outfile.write((char *)&m_altFreq[0], m_altFreq.size() * sizeof(double));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  m_outfile.write((char *)&m_maf[0], m_maf.size() * sizeof(double));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  m_outfile.write((char *)&m_avgCall[0], m_avgCall.size() * sizeof(double));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
  m_outfile.write((char *)&m_rSquared[0], m_rSquared.size() * sizeof(double));
//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
}

// Need to change the first two elements to std::vector<string> to allow for mulitple input files
int CWriteBinaryDosageFromVCF::WriteBinaryDosageFile(const std::string &vcfFilename, const std::string &infoFilename, const std::string &outputFilename) {
  m_numGroups = 1;
  // Read the information from the info file
  if (ReadInfoFile(infoFilename))
    return 1;
  // Close and clear m_infile
  m_infile.close();
  m_infile.clear();
  // Open the vcf file
  m_infile.open(vcfFilename.c_str());
  if (!m_infile.good()) {
    Rcpp::Rcerr << "Unable to open vcf file:\t" << vcfFilename << std::endl;
    return 1;
  }
  // Read in the subject IDs
  if (ReadSubjectInfo())
    return 1;
  // With subject and SNP information read it is now time to write the header
  // Open the output file
  m_outfile.open(outputFilename.c_str(), std::ios::out | std::ios::binary);
  if (!m_outfile.good()) {
    Rcpp::Rcerr << "Unable to open output file:\t" << outputFilename << std::endl;
    return 1;
  }
  // Write the header - stats about the file
  WriteHeader();
  // Write the subject data - family and subject IDs
  WriteSubjects();
  // Write the information about the SNPs
  WriteSNPInfo();
  // Read and write the dosages
  return WriteDosages();
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
int VCF2BD(const std::string &vcfFilename, const std::string &infoFilename, const std::string &outputFilename) {
  CWriteBinaryDosageFromVCF wvcf;

  return wvcf.WriteBinaryDosageFile(vcfFilename, infoFilename, outputFilename);
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
