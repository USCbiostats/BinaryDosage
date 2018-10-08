#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "BDoseWriter.h"

enum class Header4pos {
  header = 0,
    version = 4,
    numSub = 8,
    numSNPs = 12,
    numGroups = 16,
    subOptions = 20,
    snpOptions = 24,
    startSub = 28,
    startSNP = 32,
    startDosage = 36,
    startGroups = 40
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                CBDoseWriter
//
/////////////////////////////////////////////////////////////////////////////////////////////////
// Base class for writing binary dosage file

// Constructor - check if file can be opened
CBDoseWriter::CBDoseWriter(const std::string &_filename, const int _format, const int _version) : m_format(_format), m_version(_version) {
  std::ofstream outfile;
  const char blankHeader[4] = { ' ', ' ', ' ', ' '};
  const char formatHeader[4][2] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
  const char versionHeader[2][2] = {{0, 1}, {0, 2}};

  m_filename = _filename;

  m_subjectDataWritten = false;
  m_numSubjects = 0;

  m_SNPDataWritten = false;
  m_numSNPs = 0;
  m_snpsWritten = 0;

  m_startDosageData = 0;

  m_geneticDataWriter = NULL;

  if (m_format < 1 || m_format > 4 || m_version < 1 || m_version > 2) {
    m_good = false;
    return;
  }

  // Having issue opening file below - it needs to exist - this creates the file
  outfile.open(m_filename.c_str());
  if (!outfile.good()) {
    outfile.close();
    m_good = false;
    return;
  }
  outfile.close();

  m_outfile.open(m_filename.c_str(), std::ios_base::out | std::ios_base::in | std::ios_base::binary | std::ios_base::trunc);
  if (!m_outfile.good()) {
    m_good = false;
    return;
  }
  m_outfile.write(blankHeader, 4);
  m_outfile.write(formatHeader[m_format], 2);
  m_outfile.write(versionHeader[m_version], 2);
  m_good = true;
}

// Destructor - deletes geneticDataWriter and closes file
CBDoseWriter::~CBDoseWriter() {
  if (m_geneticDataWriter)
    delete m_geneticDataWriter;
  m_outfile.close();
}
// Write the subject data - check if valid
int CBDoseWriter::WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID) {
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2;

  if (m_subjectDataWritten)
    m_good = false;

  if (_SID.size() == 0)
    m_good = false;
  if (_FID.size() != 0 && _FID.size() != _SID.size()) {}
  m_good = false;

  if (!m_good)
    return 1;
  return 0;
}
// Write the SNP data - checks if data is valid
int CBDoseWriter::WriteSNPData(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                                const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                                const std::vector<std::string> &_altAllele, int _snpOptions) {
  unsigned int numSNPs;
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2, strIt3, strIt4;
  std::vector<int>::const_iterator dblIt;

  if (!m_subjectDataWritten || m_SNPDataWritten)
    m_good = false;

  numSNPs = _chromosome.size();
  if (numSNPs == 0)
    m_good = false;
  if (_snpID.size() != 0 && _snpID.size() != numSNPs)
    m_good = false;
  if (_location.size() != numSNPs || _refAllele.size() != numSNPs || _altAllele.size() != numSNPs)
    m_good = false;

  if (!m_good)
    return 1;
  return 0;
}
// Check if GeneticDataWriter can be opened
int CBDoseWriter::OpenGeneticDataWriter() {
  if (!m_subjectDataWritten || m_SNPDataWritten)
    m_good = false;
  if (m_geneticDataWriter != NULL)
    m_good = false;
  if (!m_good)
    return 1;
  return 0;
}

// Write the genetic data by call ingthe geneticDataWriter
int CBDoseWriter::WriteGeneticData(const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
  if (!m_good)
    return 1;
  if (m_geneticDataWriter == NULL) {
    m_good = false;
    return 1;
  }
  if (m_geneticDataWriter->WriteData(m_outfile, _dosage, _p0, _p1, _p2)) {
    m_good = false;
    return 1;
  }
  return 0;
}

int CBDoseWriter::Finalize() {
  const char header[4] = { 'b', 'o', 's', 'e'};

  if (!m_good || !m_subjectDataWritten || !m_SNPDataWritten || m_snpsWritten == 0 || m_snpsWritten != m_numSNPs)
    return 1;
  // Write the identifier at the beginning of the file
  m_outfile.seekp(0);
  m_outfile.write(header, 4);
  // Close the file and set good to false to prevent any further writing
  m_outfile.close();
  m_good = false;

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                CBDoseWriter1
//
/////////////////////////////////////////////////////////////////////////////////////////////////
// Class for writing binary dosage file of formats 1, 2, or 3
CBDoseWriter1::CBDoseWriter1(const std::string &_filename, const std::string &_famFilename, const std::string &_mapFilename,
                             const int _format, const int _version) : CBDoseWriter(_filename, _format, _version) {
  m_famFilename = _famFilename;
  m_mapFilename = _mapFilename;
  if (m_format < 1 || m_format > 3)
    m_good = false;
  if (m_version < 1 || m_version > 2)
    m_good = false;
}

// Create the appropriate GeneticDataWriter - If already opens returns error
int CBDoseWriter1::OpenGeneticDataWriter() {
  short scale;

  if (CBDoseWriter::OpenGeneticDataWriter())
    return 1;
  // Not checking if format and version numbers are valid. That has already been done.
  if (m_format == 1) {
    if (m_version == 1)
      m_geneticDataWriter = new CDosageDataWriter(0x7ffe, m_numSubjects);
    else
      m_geneticDataWriter = new CGeneticDataWriter1(0xfffe, m_numSubjects);
  } else if (m_format == 2) {
    if (m_version == 1)
      m_geneticDataWriter = new CDosageDataWriter(10000, m_numSubjects);
    else
      m_geneticDataWriter = new CGeneticDataWriter1(10000, m_numSubjects);
  } else {
    if (m_version == 1)
      m_geneticDataWriter = new CDosageDataWriter(10000, m_numSubjects);
    else
      m_geneticDataWriter = new CGeneticDataWriter3(10000, m_numSubjects);
  }
  return 0;
}
// Write the subject data to a text file
int CBDoseWriter1::WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID) {
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2;

  if (CBDoseWriter::WriteSubjectData(_FID, _SID))
    return 1;

  outfile.open(m_famFilename.c_str());
  if (!outfile.good()) {
    return m_good = false;
    return 1;
  }

  if (_FID.size() == 0) {
    for (strIt1 = _SID.begin(); strIt1 != _SID.end(); ++strIt1)
      outfile << *strIt1 << "\t0\t0\t9\t9" << std::endl;
  }
  else {
    strIt2 = _FID.begin();
    for (strIt1 = _SID.begin(); strIt1 != _SID.end(); ++strIt1, ++strIt2)
      outfile << *strIt1 << '\t' << *strIt2 << "\t0\t0\t9\t9" << std::endl;
  }

  if (!outfile.good()) {
    m_good = false;
  } else {
    m_subjectDataWritten = true;
    m_numSubjects = _SID.size();
    if (m_format == 3)
      m_outfile.write((char *)&m_numSubjects, sizeof(int));
    if (m_outfile.good())
      m_startDosageData = m_outfile.tellp();
    else
      m_good = false;
  }

  outfile.close();

  if (!m_good)
    return 1;
  return 0;
}
// Write the SNP data to the map file
int CBDoseWriter1::WriteSNPData(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                               const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                               const std::vector<std::string> &_altAllele, int _snpOptions) {
  unsigned int numSNPs;
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2, strIt3, strIt4;
  std::vector<int>::const_iterator dblIt;

  if (CBDoseWriter::WriteSNPData(_chromosome, _snpID, _location, _refAllele, _altAllele, _snpOptions))
    return 1;

  outfile.open(m_mapFilename.c_str());
  if (!outfile.good()) {
    return m_good = false;
    return 1;
  }

  if (_snpID.size() == 0) {
    strIt3 = _refAllele.begin();
    strIt4 = _refAllele.begin();
    dblIt = _location.begin();
    for (strIt1 = _chromosome.begin(); strIt1 != _chromosome.end(); ++strIt1, ++strIt3, ++strIt4, ++dblIt)
      outfile << *strIt1 << '\t' << *strIt1 << ':' << *dblIt << "\t0\t" << *dblIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
  }
  else {
    strIt2 = _snpID.begin();
    strIt3 = _refAllele.begin();
    strIt4 = _refAllele.begin();
    dblIt = _location.begin();
    for (strIt1 = _chromosome.begin(); strIt1 != _chromosome.end(); ++strIt1, ++strIt2, ++strIt3, ++strIt4, ++dblIt)
      outfile << *strIt1 << '\t' << *strIt2 << "\t0\t" << *dblIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
  }

  if (!outfile.good()) {
    m_good = false;
  } else {
    m_SNPDataWritten = true;
    m_numSNPs = _chromosome.size();
  }

  outfile.close();

  if (!m_good)
    return 1;

  return OpenGeneticDataWriter();
}
//////////////////////////////////////////////////////////////////////////////////////////////////
//
//                                CBDoseWriter4
//
/////////////////////////////////////////////////////////////////////////////////////////////////
// Class for writing binary dosage file of format 4
CBDoseWriter4::CBDoseWriter4(const std::string &_filename, const int _format,
                             const int _version) : CBDoseWriter(_filename, _format, _version) {
  if (m_format == 4 || m_version < 1 || m_version > 2)
    m_good = false;
}
// Process a vector of strings and write to the file
int CBDoseWriter4::ProcessString(const std::vector<std::string> &_stringsToWrite, unsigned int &_strLength) {
  std::vector<std::string>::const_iterator strIt;
  std::ostringstream oss;
  std::string s1;

  oss.str("");
  for (strIt = _stringsToWrite.begin(); strIt != _stringsToWrite.end(); ++strIt)
    oss << *strIt << '\t';
  s1 = oss.str();
  s1[s1.length() - 1] = 0;
  m_outfile.write(s1.c_str(), s1.length());
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }
  _strLength = s1.length();

  return 0;
}

// Process a vector of double vectors and write to file
int CBDoseWriter4::ProcessExtraSNPData(const std::vector<std::vector<double> > &_dataToWrite, const unsigned int _groupSize) {
  std::vector<std::vector<double> >::const_iterator vdblIt;

  for (vdblIt = _dataToWrite.begin(); vdblIt != _dataToWrite.end(); ++vdblIt) {
    if (vdblIt->size() != _groupSize)
      return 1;
    m_outfile.write((char *)vdblIt->data(), _groupSize * sizeof(double));
    if (m_outfile.good()) {
      m_good = false;
      return 1;
    }
  }
  return 0;
}
// Figure out what SNP options are used
unsigned int CBDoseWriter4::GetSNPOptions(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                                          const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                                          const std::vector<std::string> &_altAllele,
                                          const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
                                          const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq) {
  unsigned int retVal = 0x0014;
  std::string s1;
  std::vector<std::string>::const_iterator strIt;

  if (_snpID.size() != 0 && _snpID[0] != "")
    retVal |= 0x0002;
  s1 = _chromosome[0];
  for (strIt = _chromosome.begin(); strIt != _chromosome.end(); ++strIt) {
    if (*strIt != s1)
      break;
  }
  if (strIt == _chromosome.end())
    retVal |= 0x0008;
  if (_refAllele.size() != 0)
    retVal |= 0x0020;
  if (_altAllele.size() != 0)
    retVal |= 0x0040;
  if (_altFreq.size() != 0)
    retVal |= 0x0080;
  if (_maf.size() != 0)
    retVal |= 0x0100;
  if (_avgCall.size() != 0)
    retVal |= 0x0200;
  if (_rSq.size() != 0)
    retVal |= 0x0400;
  return retVal;
}
