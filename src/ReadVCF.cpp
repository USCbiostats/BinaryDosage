#include <string>
#include <fstream>
#include <sstream>
#include <Rcpp.h>
#include "ReadVCF.h"

const double NumberValue[4][10] = {
  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
  {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9},
  {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09},
  {0.000, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009}
};

CReadVCFBase::CReadVCFBase(const std::string &filename) {
  m_filename = filename;
  m_infile.open(m_filename.c_str());
  
  m_subjectLine = 0;
  m_startLine = 0;
  
  m_numSubjects = 0;
  m_numSNPs = 0;
}

CReadVCFBase::~CReadVCFBase() {
  m_infile.close();
}

int CReadVCFBase::ReadSubjects() {
  const std::string columnNames[9] = { "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"};
  std::string junk;
  std::string colName;
  std::string subjectID;
  std::istringstream iss;
  int i;
  
  if (!m_infile.good())
    return 1;

  std::getline(m_infile, junk);
  if (junk.substr(0, 16) != "##fileformat=VCF")
    return 2;
  
  m_subjectLine = 1;
  do {
    std::getline(m_infile, junk);
    ++m_subjectLine;
  } while (m_infile.good() && junk.substr(0,2) == "##");
  
  if (!m_infile.good())
    return 3;
  if (junk.substr(0,6) != "#CHROM")
    return 4;
  
  iss.str(junk);
  for (i = 0; i < 9; ++i) {
    iss >> colName;
    if (colName != columnNames[i])
      return 5;
  }
  
  m_numSubjects = 0;
  do {
    iss >> subjectID;
    if (iss.fail())
      break;
    ++m_numSubjects;
  } while (1);
  
  iss.clear();
  iss.seekg(0);
  for (i = 0; i < 9; ++i)
    iss >> colName;
  m_subjectID.resize(m_numSubjects);
  for (std::vector<std::string>::iterator sub = m_subjectID.begin(); sub != m_subjectID.end(); ++sub) {
    iss >> subjectID;
    *sub = subjectID;
  }
  
  m_startLine = m_subjectLine + 1;
  
  return 0;
}

double CReadVCFBase::AlternateAlleleFrequency() {
  if (m_dosage.size() == 0)
    return 0.;
  return std::accumulate(m_dosage.begin(), m_dosage.end(), 0.) / (2. * m_dosage.size());
}
// ******************************************************************************************
//                                CReadVCF_HRC
// ******************************************************************************************

double CReadVCF_HRC::ReadDosage(const char *x) {
  double sum;
  
  if (*x < '0' || *x > '2')
    return -1;
  sum = NumberValue[0][*x - '0'];
  
  ++x;
  if (*x != '.')
    return -1;
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[1][*x - '0'];
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[2][*x - '0'];
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[3][*x - '0'];
  
  return sum;
}

double CReadVCF_HRC::ReadProbability(const char *x) {
  double sum;
  
  if (*x < '0' || *x > '1')
    return -1;
  sum = NumberValue[0][*x - '0'];
  
  ++x;
  if (*x != '.')
    return -1;
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[1][*x - '0'];
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[2][*x - '0'];
  
  ++x;
  if (*x < '0' || *x > '9')
    return -1;
  sum += NumberValue[3][*x - '0'];
  
  return sum;
}

int CReadVCF_HRC::GetFirstSNP() {
  std::string junk;
  std::string skipped;
  std::string valueString;
  std::istringstream iss;
  std::vector<double>::iterator d, p0, p1, p2;
  int i;

  m_numSNPs = 0;  
  
  if (!m_infile.is_open())
    return 1;

  if (m_numSubjects == 0)
    return 1;
  
  m_snpID.resize(1);
  m_chromosome.resize(1);
  m_bp.resize(1);
  m_refAllele.resize(1);
  m_altAllele.resize(1);
  m_dosage.resize(m_numSubjects);
  m_p0.resize(m_numSubjects);
  m_p1.resize(m_numSubjects);
  m_p2.resize(m_numSubjects);

  m_infile.clear();
  m_infile.seekg(0);
  for (i = 1; i < m_startLine; ++i)
    getline(m_infile, junk);

  return GetNextSNP();  
}

int CReadVCF_HRC::GetNextSNP() {
  std::string junk;
  std::string skipped;
  std::string valueString;
  std::istringstream iss;
  std::vector<double>::iterator d, p0, p1, p2;
  int i;
  
  if (!m_infile.is_open())
    return 1;
  
  if (!m_infile.good())
    return 1;
  
  if (m_numSubjects == 0)
    return 1;
  
  if (m_dosage.size() == 0)
    return 1;

  getline(m_infile, junk);
  if (m_infile.fail())
    return 1;
  
  ++m_numSNPs;
  iss.str(junk);
  iss >> m_chromosome[0] >> m_bp[0] >> m_snpID[0] >> m_refAllele[0] >> m_altAllele[0];
  for (i = 0; i < 4; ++i)
    iss >> skipped;
  
  for (d = m_dosage.begin(), p0 = m_p0.begin(), p1 = m_p1.begin(), p2 = m_p2.begin(); d != m_dosage.end(); ++d, ++p0, ++p1, ++p2) {
    iss >> valueString;
    *d = ReadDosage(&valueString[4]);
    *p0 = ReadProbability(&valueString[10]);
    *p1 = ReadProbability(&valueString[16]);
    *p2 = ReadProbability(&valueString[22]);
  }
  
  return 0;
}