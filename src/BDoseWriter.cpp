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

  m_groupDataWritten = false;
  m_numGroups = 0;

  m_subjectDataWritten = false;
  m_numSubjects = 0;

  m_SNPDataWritten = false;
  m_numSNPs = 0;
  m_snpsWritten = 0;

  m_startDosageData = 0;

  m_geneticDataWriter = NULL;

  // Is the format valid
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
  m_outfile.write(formatHeader[m_format - 1], 2);
  m_outfile.write(versionHeader[m_version - 1], 2);
  m_good = true;
}

// Destructor - deletes geneticDataWriter and closes file
CBDoseWriter::~CBDoseWriter() {
  if (m_geneticDataWriter)
    delete m_geneticDataWriter;
  m_outfile.close();
}
// Check if GeneticDataWriter can be opened - Have required steps been done - Opening is done in derived class
int CBDoseWriter::OpenGeneticDataWriter() {
  if (!m_subjectDataWritten || !m_SNPDataWritten)
    m_good = false;
  if (m_geneticDataWriter != NULL)
    m_good = false;
  if (!m_good)
    return 1;
  return 0;
}
// Check if group data is valid
int CBDoseWriter::WriteGroupData(const std::vector<int> &_groupSizes) {
  m_numGroups = _groupSizes.size();
  if (m_numGroups == 0) {
    m_good = false;
    return 1;
  }
  m_groupDataWritten = true;
  return 0;
}
// Write the subject data - check if valid - writing is done in the derived class
int CBDoseWriter::WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID) {
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2;

  if (!m_groupDataWritten || m_subjectDataWritten)
    m_good = false;

  if (_SID.size() == 0)
    m_good = false;
  if (_FID.size() != 0 && _FID.size() != _SID.size())
    m_good = false;

  if (!m_good)
    return 1;

  m_numSubjects = _SID.size();
  return 0;
}
// Write the SNP data - checks if data is valid - Writing is done in the derived class
int CBDoseWriter::WriteSNPData(const std::vector<std::string> &_chromosome,
                               const std::vector<std::string> &_snpID,
                               const std::vector<int> &_location,
                               const std::vector<std::string> &_refAllele,
                               const std::vector<std::string> &_altAllele,
                               const std::vector<std::vector<double> > &_aaf,
                               const std::vector<std::vector<double> > &_maf,
                               const std::vector<std::vector<double> > &_avgCall,
                               const std::vector<std::vector<double> > &_rSq) {
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2, strIt3, strIt4;
  std::vector<int>::const_iterator dblIt;
  unsigned int numSNPs;

  if (!m_groupDataWritten || !m_subjectDataWritten || m_SNPDataWritten)
    m_good = false;

  numSNPs = _location.size();
  m_numSNPs = numSNPs;
  if (m_numSNPs == 0)
    m_good = false;
  if (_chromosome.size() != 1 && _chromosome.size() != numSNPs)
    m_good = false;
  if (_snpID.size() != 0 && _snpID.size() != numSNPs)
    m_good = false;
  if (_refAllele.size() != numSNPs || _altAllele.size() != numSNPs)
    m_good = false;

  if (!m_good)
    return 1;
  return 0;
}

// Write the genetic data by calling the geneticDataWriter - Should be same for all derived classes
int CBDoseWriter::WriteGeneticData(const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
  if (!m_good)
    return 1;
  if (!m_groupDataWritten || !m_subjectDataWritten || !m_SNPDataWritten) {
    m_good = false;
    return 1;
  }
  if (m_geneticDataWriter == NULL) {
    m_good = false;
    return 1;
  }
  if (m_geneticDataWriter->WriteData(m_outfile, _dosage, _p0, _p1, _p2)) {
    m_good = false;
    return 1;
  }
  ++m_snpsWritten;
  return 0;
}
// Test if required steps have been done - Should be same for all classes
int CBDoseWriter::Finalize() {
  const char header[4] = { 'b', 'o', 's', 'e'};

  if (!m_good || !m_groupDataWritten || !m_subjectDataWritten || !m_SNPDataWritten || m_snpsWritten == 0 || m_snpsWritten != m_numSNPs)
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
  // Check if version is valid
  if (m_format < 1 || m_format > 3)
    m_good = false;
  if (m_version < 1 || m_version > 2)
    m_good = false;
}

// Create the appropriate GeneticDataWriter - If already opens returns error
int CBDoseWriter1::OpenGeneticDataWriter() {
  if (CBDoseWriter::OpenGeneticDataWriter())
    return 1;
  // Check if format and version are good. Shouldn't be required but who knows?
  if (m_format == 1) {
    if (m_version == 1) {
      m_geneticDataWriter = new CDosageDataWriter(0x7ffe, m_numSubjects);
    } else if (m_version == 2) {
      m_geneticDataWriter = new CGeneticDataWriter1(0xfffe, m_numSubjects);
    } else {
      m_good = false;
      return 1;
    }
  } else if (m_format == 2) {
    if (m_version == 1)
      m_geneticDataWriter = new CDosageDataWriter(10000, m_numSubjects);
    else
      m_geneticDataWriter = new CGeneticDataWriter1(10000, m_numSubjects);
  } else if (m_format == 3) {
    if (m_version == 1) {
      m_geneticDataWriter = new CDosageDataWriter(10000, m_numSubjects);
    } else if (m_version == 2) {
      m_geneticDataWriter = new CGeneticDataWriter3(10000, m_numSubjects);
    } else {
      m_good = false;
      return 1;
    }
  } else {
    m_good = false;
    return 1;
  }
  return 0;
}
// Write the subject data to a text file - This is a fam file in plink format - Slightly modified allow subject IDs only
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

  if (_FID.size() == 0 || _FID[0] == "") {
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
    m_numSubjects = _SID.size();
    if (m_format == 3)
      m_outfile.write((char *)&m_numSubjects, sizeof(int));
    if (m_outfile.good()) {
      m_startDosageData = m_outfile.tellp();
      m_subjectDataWritten = true;
    } else {
      m_good = false;
    }
  }

  outfile.close();

  if (!m_good)
    return 1;
  return 0;
}
// Write the SNP data to the map file - This is in an extened map file format used by plink
int CBDoseWriter1::WriteSNPData(const std::vector<std::string> &_chromosome,
                                const std::vector<std::string> &_snpID,
                                const std::vector<int> &_location,
                                const std::vector<std::string> &_refAllele,
                                const std::vector<std::string> &_altAllele,
                                const std::vector<std::vector<double> > &_aaf,
                                const std::vector<std::vector<double> > &_maf,
                                const std::vector<std::vector<double> > &_avgCall,
                                const std::vector<std::vector<double> > &_rSq) {
  std::ofstream outfile;
  std::vector<std::string>::const_iterator strIt1, strIt2, strIt3, strIt4;
  std::vector<int>::const_iterator intIt;

  if (CBDoseWriter::WriteSNPData(_chromosome, _snpID, _location, _refAllele, _altAllele, _aaf, _maf, _avgCall, _rSq))
    return 1;
  outfile.open(m_mapFilename.c_str());
  if (!outfile.good()) {
    return m_good = false;
    return 1;
  }
  if (_snpID.size() == 0) {
    strIt1 = _chromosome.begin();
    strIt3 = _refAllele.begin();
    strIt4 = _refAllele.begin();
    for (intIt = _location.begin(); intIt != _location.end(); ++strIt3, ++strIt4, ++intIt) {
      outfile << *strIt1 << '\t' << *strIt1 << ':' << *intIt << "\t0\t" << *intIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
      if (_chromosome.size() != 1)
        ++strIt1;
    }
  }
  else {
    strIt1 = _chromosome.begin();
    strIt2 = _snpID.begin();
    strIt3 = _refAllele.begin();
    strIt4 = _refAllele.begin();
    for (intIt = _location.begin(); intIt != _location.end(); ++strIt2, ++strIt3, ++strIt4, ++intIt) {
      outfile << *strIt1 << '\t' << *strIt2 << "\t0\t" << *intIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
      if (_chromosome.size() != 1)
        ++strIt1;
    }
  }

  if (!outfile.good())
    m_good = false;
  else
    m_SNPDataWritten = true;

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
  const int blankValues[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

  m_startSubjectData = 0;
  m_startSNPData = 0;
  m_startSNPInfo = 0;
  m_startDosageData = 0;
  m_snpOptions = 0;

  if (m_format != 4 || m_version < 1 || m_version > 2) {
    m_good = false;
    return;
  }
  // Write out a blank header
  m_outfile.write((const char *)blankValues, 8 * sizeof(int));
  if (!m_outfile.good())
    m_good = false;
}

// Process a vector of strings and write to the file
int CBDoseWriter4::ProcessString(const std::vector<std::string> &_stringsToWrite, int &_strLength) {
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
int CBDoseWriter4::ProcessExtraSNPData(const std::vector<std::vector<double> > &_dataToWrite, const int _groupSize) {
  std::vector<std::vector<double> >::const_iterator vdblIt;

  for (vdblIt = _dataToWrite.begin(); vdblIt != _dataToWrite.end(); ++vdblIt) {
    if ((int)vdblIt->size() != _groupSize) {
      m_good = false;
      return 1;
    }
    m_outfile.write((char *)vdblIt->data(), _groupSize * sizeof(double));
    if (!m_outfile.good()) {
      m_good = false;
      return 1;
    }
  }
  return 0;
}
// Figure out what SNP options are used
// Unsure if this should be done or if user should provide the options and then have the
// program check if the options are valid ???
int CBDoseWriter4::GetSNPOptions(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                                 const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                                 const std::vector<std::string> &_altAllele,
                                 const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
                                 const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq) {
  unsigned int retVal = 0x0014;

  if (_snpID.size() != 0)
    retVal |= 0x0002;

  if (_chromosome.size() == 1)
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
// Create the appropriate GeneticDataWriter - If already opens returns error
int CBDoseWriter4::OpenGeneticDataWriter() {
  if (CBDoseWriter::OpenGeneticDataWriter())
    return 1;
  // Check if format and version are good. Shouldn't be required but who knows?
  if (m_format == 4) {
    if (m_version == 1) {
      m_geneticDataWriter = new CDosageDataWriter(10000, m_numSubjects);
    } else if (m_version == 2) {
      m_geneticDataWriter = new CGeneticDataWriter3(10000, m_numSubjects);
    } else {
      m_good = false;
      return 1;
    }
  } else {
    m_good = false;
    return 1;
  }
  return 0;
}

// Write out the group data
int CBDoseWriter4::WriteGroupData(const std::vector<int> &_groupSizes) {
  int startSubjects;
  if (CBDoseWriter::WriteGroupData(_groupSizes))
    return 1;

  m_groupDataWritten = false;
  m_outfile.seekp((int)Header4pos::numGroups);
  m_outfile.write((char *)&m_numGroups, sizeof(int));
  m_outfile.seekp((int)Header4pos::startGroups);
  m_outfile.write((char *)_groupSizes.data(), m_numGroups * sizeof(int));
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }
  m_startSubjectData = m_outfile.tellp();
  startSubjects = (int) m_startSubjectData;
  m_outfile.seekp((int)Header4pos::startSub);
  m_outfile.write((char *)&startSubjects, sizeof(int));
  m_outfile.seekp(m_startSubjectData);
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  m_groupDataWritten = true;
  return 0;
}
// Write the subject data
int CBDoseWriter4::WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID) {
  const int zeroInt = 0;
  int sidLength, fidLength;
  int startSNP;

  if (CBDoseWriter::WriteSubjectData(_FID, _SID))
    return 1;

  m_numSubjects = _SID.size();
  sidLength = 0;
  fidLength = 0;

  m_outfile.seekp(m_startSubjectData);
  m_outfile.write((char *)&zeroInt, sizeof(int));
  m_outfile.write((char *)&zeroInt, sizeof(int));
  if (ProcessString(_SID, sidLength)) {
    m_good = false;
    return 1;
  }
  if (_FID.size() == 0 || _FID[0] == "") {
    fidLength = 0;
  } else {
    if (ProcessString(_FID, fidLength)) {
      m_good = false;
      return 1;
    }
  }

  m_startSNPData = m_outfile.tellp();
  startSNP = (int)m_startSNPData;

  m_outfile.seekp((int)Header4pos::numSub);
  m_outfile.write((char *)&m_numSubjects, sizeof(int));
  m_outfile.seekp(m_startSubjectData);
  m_outfile.write((char *)&sidLength, sizeof(int));
  m_outfile.write((char *)&fidLength, sizeof(int));
  m_outfile.seekp((int)Header4pos::startSNP);
  m_outfile.write((char *)&startSNP, sizeof(int));
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  m_outfile.seekp(m_startSNPData);
  m_subjectDataWritten = true;
  return 0;
}

int CBDoseWriter4::WriteSNPData(const std::vector<std::string> &_chromosome,
                                const std::vector<std::string> &_snpID,
                                const std::vector<int> &_location,
                                const std::vector<std::string> &_refAllele,
                                const std::vector<std::string> &_altAllele,
                                const std::vector<std::vector<double> > &_aaf,
                                const std::vector<std::vector<double> > &_maf,
                                const std::vector<std::vector<double> > &_avgCall,
                                const std::vector<std::vector<double> > &_rSq) {
  const int zeroInt = 0;
  const char zeroChar = 0;
  int snpLength, chrLength, refLength, altLength;
  int snpOptions;
  int startDosage;
  std::string s1;

  if (CBDoseWriter::WriteSNPData(_chromosome, _snpID, _location, _refAllele, _altAllele, _aaf, _maf, _avgCall, _rSq))
    return 1;
  snpOptions = GetSNPOptions(_chromosome, _snpID,  _location, _refAllele, _altAllele, _aaf, _maf, _avgCall, _rSq);
  m_snpOptions = snpOptions;
  m_outfile.seekp((int)Header4pos::snpOptions);
  m_outfile.write((char *)&snpOptions, sizeof(int));
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  m_outfile.seekp(m_startSNPData);
  m_outfile.write((char *)&zeroInt, sizeof(int));
  m_outfile.write((char *)&zeroInt, sizeof(int));
  m_outfile.write((char *)&zeroInt, sizeof(int));
  m_outfile.write((char *)&zeroInt, sizeof(int));
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  if (_snpID.size() == 0) {
    snpLength = 0;
  } else {
    if (ProcessString(_snpID, snpLength))
      return 1;
  }

  if (snpOptions & 0x0008) {
    s1 = _chromosome[0];
    m_outfile.write(s1.c_str(), s1.length());
    m_outfile.write(&zeroChar, 1);
    chrLength = s1.length() + 1;
    if (!m_outfile.good()) {
      m_good = false;
      return 1;
    }
  }
  else {
    if (ProcessString(_chromosome, chrLength))
      return 1;
  }

  m_outfile.write((char *)_location.data(), m_numSNPs * sizeof(unsigned int));
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  if (_refAllele.size() != 0) {
    if (ProcessString(_refAllele, refLength))
      return 1;
  } else {
    refLength = 0;
  }
  if (_altAllele.size() != 0) {
    if (ProcessString(_altAllele, altLength))
      return 1;
  } else {
    altLength = 0;
  }

  m_startSNPInfo = m_outfile.tellp();
  if (snpOptions & 0x0080) {
    if (ProcessExtraSNPData(_aaf, m_numGroups))
      return 1;
  }
  if (snpOptions & 0x0100) {
    if (ProcessExtraSNPData(_maf, m_numGroups))
      return 1;
  }
  if (snpOptions & 0x0200) {
    if (ProcessExtraSNPData(_avgCall, m_numGroups))
      return 1;
  }
  if (snpOptions & 0x0400) {
    if (ProcessExtraSNPData(_rSq, m_numGroups))
      return 1;
  }

  m_startDosageData = m_outfile.tellp();
  startDosage = (int)m_startDosageData;

  m_outfile.seekp((int)Header4pos::numSNPs);
  m_outfile.write((char *)&m_numSNPs, sizeof(int));
  m_outfile.seekp(m_startSNPData);
  m_outfile.write((char *)&snpLength, sizeof(int));
  m_outfile.write((char *)&chrLength, sizeof(int));
  m_outfile.write((char *)&refLength, sizeof(int));
  m_outfile.write((char *)&altLength, sizeof(int));

  m_outfile.seekp((int)Header4pos::startDosage);
  m_outfile.write((char *)&startDosage, sizeof(int));
  m_outfile.seekp(m_startDosageData);

  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  m_SNPDataWritten = true;
  return OpenGeneticDataWriter();
}

int CBDoseWriter4::UpdateSNPInfo(const std::vector<std::vector<double> > &_aaf,
                                 const std::vector<std::vector<double> > &_maf,
                                 const std::vector<std::vector<double> > &_avgCall,
                                 const std::vector<std::vector<double> > &_rSq) {
  if (!m_good)
    return 1;

  if (!m_groupDataWritten || !m_subjectDataWritten || !m_SNPDataWritten) {
    m_good = false;
    return 1;
  }

  m_outfile.seekp(m_startSNPInfo);
  if (m_snpOptions & 0x0080) {
    if (ProcessExtraSNPData(_aaf, m_numGroups))
      return 1;
  }
  if (m_snpOptions & 0x0100) {
    if (ProcessExtraSNPData(_maf, m_numGroups))
      return 1;
  }
  if (m_snpOptions & 0x0200) {
    if (ProcessExtraSNPData(_avgCall, m_numGroups))
      return 1;
  }
  if (m_snpOptions & 0x0400) {
    if (ProcessExtraSNPData(_rSq, m_numGroups))
      return 1;
  }
  if (!m_outfile.good()) {
    m_good = false;
    return 1;
  }

  return 0;
}
