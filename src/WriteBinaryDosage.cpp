#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "WriteBinaryDosage.h"
#include <Rcpp.h>

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
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage::CWriteBinaryDosage(const std::vector<std::string> &filenames) {
  memset(m_version, 0, sizeof(m_version));
  if (filenames.size() == 0)
    m_ready = false;
  m_filename = filenames[0];
  if (m_filename != "") {
    m_outfile.open(m_filename.c_str(), std::ios_base::out | std::ios_base::binary);
    if (m_outfile.good())
      m_ready = true;
  }
}

CWriteBinaryDosage::~CWriteBinaryDosage() {
  m_outfile.close();
}

short CWriteBinaryDosage::ConvertToShort(const double x, const double scale) {
  short s1, s2;

  s1 = x * scale;
  s2 = s1 + 1;
  if (fabs(x - (double)s1 / scale) < fabs(x - (double)s2 / scale))
    return s1;
  return s2;
}

int CWriteBinaryDosage::WriteHeader() {
  const char header[4] = { 'b', 'o', 's', 'e' };
  if (!m_ready)
    return 1;

  m_outfile.write(header, 4);
  m_outfile.write(m_version, 4);
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage::AddDosagesOnly(const std::vector<std::vector<double> > &dosageValues, double scale) {
  std::vector<double>::const_iterator doseIt;
  std::vector<short>::iterator sdIt;

  if (!m_ready)
    return 1;
  if (dosageValues.size() != 1 && dosageValues.size() != 4) {
    m_ready = false;
    return 1;
  }
  if (dosageValues[0].size() == 0) {
    m_ready = false;
    return 1;
  }
  if (m_dataToWrite.size() != dosageValues[0].size()) {
    if (m_dataToWrite.size() != 0) {
      m_ready = false;
      return 1;
    }
    m_dataToWrite.resize(dosageValues[0].size());
  }

  sdIt = m_dataToWrite.begin();
  for (doseIt = dosageValues[0].begin(); doseIt != dosageValues[0].end(); ++doseIt, ++sdIt)
    *sdIt = ConvertToShort(*doseIt, scale);
  m_outfile.seekp(0, m_outfile.end);
  m_outfile.write((char *)m_dataToWrite.data(), m_dataToWrite.size() * sizeof(short));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage::AddGeneticValues32(const std::vector<std::vector<double> > &geneticValues) {
  const double scale = 10000.;
  const short shortScale = 10000;
  short dose, p0, p1, p2;
  std::vector<double>::const_iterator doseIt, p0It, p1It, p2It;
  std::vector<short>::iterator sdIt, extraIt;
  int numExtra = 0;
  int writeSize;

  if (!m_ready) {
    Rcpp::Rcout << "Error already exists" << std::endl;
    return 1;
  }
  if (geneticValues.size() != 4) {
    Rcpp::Rcout << "Wrong size" << std::endl;
    m_ready = false;
    return 1;
  }
  if (geneticValues[0].size() == 0) {
    Rcpp::Rcout << "No values" << std::endl;
    m_ready = false;
    return 1;
  }
  if (geneticValues[0].size() != geneticValues[1].size() || geneticValues[0].size() != geneticValues[2].size() || geneticValues[0].size() != geneticValues[3].size()) {
    Rcpp::Rcout << "Different sizes" << std::endl;
    m_ready = false;
    return 1;
  }
  if (m_dataToWrite.size() != 4 * geneticValues[0].size()) {
    if (m_dataToWrite.size() != 0) {
      m_ready = false;
      return 1;
    }
    m_dataToWrite.resize(4 * geneticValues[0].size());
  }

  doseIt = geneticValues[0].begin();
  p0It = geneticValues[1].begin();
  p1It = geneticValues[2].begin();
  p2It = geneticValues[3].begin();
  sdIt = m_dataToWrite.begin();
  extraIt = sdIt + geneticValues[0].size();
  numExtra = 0;
  for (; doseIt != geneticValues[0].end(); ++doseIt, ++p0It, ++p1It, ++p2It, ++sdIt) {
    dose = ConvertToShort(*doseIt, scale);
    p0 = ConvertToShort(*p0It, scale);
    p1 = ConvertToShort(*p1It, scale);
    p2 = ConvertToShort(*p2It, scale);
    if (p0 + p1 + p2 != shortScale || p1 + p2 + p2 != dose) {
      *sdIt = dose | 0x8000;
      *extraIt = p1 | 0x8000;
      ++extraIt;
      *extraIt = p0;
      ++extraIt;
      *extraIt = p2;
      ++extraIt;
      numExtra += 3;
    } else if (p0 != 0 && p2 != 0){
      *sdIt = dose | 0x8000;
      *extraIt = p1;
      ++extraIt;
      ++numExtra;
    } else {
      *sdIt = dose;
    }
  }
  writeSize = (geneticValues[0].size() + numExtra) * sizeof(short);
  m_outfile.seekp(0, m_outfile.end);
  m_outfile.write((char *)&writeSize, sizeof(int));
  m_outfile.write((char *)m_dataToWrite.data(), writeSize);
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteMultifileBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CWriteMultifileBinaryDosage::CWriteMultifileBinaryDosage(const std::vector<std::string> &filenames) : CWriteBinaryDosage(filenames) {
  if (m_ready) {
    m_ready = false;
    if (filenames.size() == 3) {
      m_famFilename = filenames[1];
      m_mapFilename = filenames[2];
      m_famFile.open(m_famFilename.c_str());
      m_mapFile.open(m_mapFilename.c_str());
      if (m_famFile.good() && m_mapFile.good())
        m_ready = true;
    }
  }
}

int CWriteMultifileBinaryDosage::AddGeneticValues1or2(const std::vector<std::vector<double> > &geneticValues, double scale) {
  std::vector<double>::const_iterator gp0It, gp1It;
  std::vector<short>::iterator sp0It, sp1It;

  if (!m_ready)
    return 1;
  if (geneticValues.size() != 4) {
    m_ready = false;
    return 1;
  }
  if (geneticValues[2].size() == 0) {
    m_ready = false;
    return 1;
  }
  if (geneticValues[2].size() != geneticValues[3].size()) {
    m_ready = false;
    return 1;
  }
  if (m_dataToWrite.size() != 2 * geneticValues[2].size()) {
    if (m_dataToWrite.size() != 0) {
      m_ready = false;
      return 1;
    }
    m_dataToWrite.resize(2 * geneticValues[2].size());
  }

  sp0It = m_dataToWrite.begin();
  sp1It = sp0It + geneticValues[2].size();
  gp0It = geneticValues[2].begin();
  gp1It = geneticValues[3].begin();
  for (; gp0It != geneticValues[2].end(); ++gp0It, ++gp1It, ++sp0It, ++sp1It) {
    *sp0It = ConvertToShort(*gp0It, scale);
    *sp1It = ConvertToShort(*gp1It, scale);
  }
  m_outfile.seekp(0, m_outfile.end);
  m_outfile.write((char *)m_dataToWrite.data(), m_dataToWrite.size() * sizeof(short));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

CWriteMultifileBinaryDosage::~CWriteMultifileBinaryDosage() {
  m_famFile.close();
  m_mapFile.close();
}

int CWriteMultifileBinaryDosage::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  std::vector<std::string>::const_iterator iid, fid;

  if (!m_ready)
    return 1;

  if (FID.size() == 0) {
    for (std::vector<std::string>::const_iterator iid = IID.begin(); iid != IID.end(); ++iid)
      m_famFile << *iid << "\t0\t0\t9\t9" << std::endl;
  } else {
    for (iid = IID.begin(), fid = FID.begin(); iid != IID.end(); ++iid, ++fid)
      m_famFile << *fid << '\t' << *iid << "\t0\t0\t9\t9" << std::endl;
  }
  if (!m_famFile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteMultifileBinaryDosage::AddSNP(const std::string &chromosome, const std::string &snpID, int location,
                                          const std::string &refAllele, const std::string &altAllele,
                                          const std::vector<double> &altFreq, const std::vector<double> &maf,
                                          const std::vector<double> &avgCall, const std::vector<double> &rSq) {
  if (!m_ready)
    return 1;
  m_mapFile << chromosome << '\t' << snpID << "\t0\t" << location << '\t' << refAllele << '\t' << altAllele << std::endl;
  if (!m_mapFile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteMultifileBinaryDosage::WriteAllSNPs(const std::vector<std::string> &chromosome,
                         const std::vector<std::string> &snpID,
                         const std::vector<int> bp,
                         const std::vector<std::string> &refAllele,
                         const std::vector<std::string> &altAllele,
                         const std::vector<std::vector<double> > &altFreq,
                         const std::vector<std::vector<double> > &maf,
                         const std::vector<std::vector<double> > &avgCall,
                         const std::vector<std::vector<double> > &rSq) {
  unsigned int expSize;
  bool generateSNPID, singleChromosome;
  std::vector<std::string>::const_iterator chrIt, snpIt, refIt, altIt;
  std::vector<int>::const_iterator bpIt;

  if (!m_ready)
    return 1;
  expSize = bp.size();
  if (expSize == 0) {
    m_ready = false;
    return 1;
  };
  generateSNPID = (snpID.size() == 0);
  if (!generateSNPID && snpID.size() != expSize) {
    m_ready = false;
    return 1;
  }
  singleChromosome = (chromosome.size() == 1);
  if (!singleChromosome && chromosome.size() != expSize) {
    m_ready = false;
    return 1;
  }
  if (refAllele.size() != expSize || altAllele.size() != 1) {
    m_ready = false;
    return 1;
  }

  chrIt = chromosome.begin();
  bpIt = bp.begin();
  refIt = refAllele.begin();
  altIt = altAllele.begin();
  if (generateSNPID) {
    for (;bpIt != bp.end(); ++bpIt, ++refIt, ++altIt) {
      m_mapFile << *chrIt << '\t' << *chrIt << ':' << *bpIt << "\t0\t" << *bpIt << '\t' << *refIt << '\t' << *altIt << std::endl;
      if (!singleChromosome)
        ++chrIt;
    }
  } else {
    snpIt = snpID.begin();
    for (;bpIt != bp.end(); ++bpIt, ++snpIt, ++refIt, ++altIt) {
      m_mapFile << *chrIt << '\t' << *snpIt << "\t0\t" << *bpIt << '\t' << *refIt << '\t' << *altIt << std::endl;
      if (!singleChromosome)
        ++chrIt;
    }
  }
  if (!m_mapFile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage11
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage11::CWriteBinaryDosage11(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x01;
  m_version[3] = 0x01;
}

int CWriteBinaryDosage11::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddDosagesOnly(geneticValues, 32767.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage12
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage12::CWriteBinaryDosage12(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x01;
  m_version[3] = 0x02;
}

int CWriteBinaryDosage12::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddGeneticValues1or2(geneticValues, 65534.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage21
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage21::CWriteBinaryDosage21(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x02;
  m_version[3] = 0x01;
}

int CWriteBinaryDosage21::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddDosagesOnly(geneticValues, 10000.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage22
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage22::CWriteBinaryDosage22(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x02;
  m_version[3] = 0x02;
}

int CWriteBinaryDosage22::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddGeneticValues1or2(geneticValues, 10000.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage31
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage31::CWriteBinaryDosage31(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x03;
  m_version[3] = 0x01;
}

int CWriteBinaryDosage31::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  if (!m_ready)
    return 1;
  numSubjects = IID.size();
  m_outfile.write((char *)&numSubjects, sizeof(int));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return CWriteMultifileBinaryDosage::WriteSubjects(FID, IID);
}

int CWriteBinaryDosage31::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddDosagesOnly(geneticValues, 10000.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage32
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage32::CWriteBinaryDosage32(const std::vector<std::string> &filenames) : CWriteMultifileBinaryDosage(filenames) {
  m_version[1] = 0x03;
  m_version[3] = 0x02;
}

int CWriteBinaryDosage32::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  if (!m_ready)
    return 1;
  numSubjects = IID.size();
  m_outfile.write((char *)&numSubjects, sizeof(int));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return CWriteMultifileBinaryDosage::WriteSubjects(FID, IID);
}

int CWriteBinaryDosage32::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddGeneticValues32(geneticValues);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage4x
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage4x::CWriteBinaryDosage4x(const std::vector<std::string> &filenames) : CWriteBinaryDosage(filenames) {
  std::string outputFilename;
  if (m_ready) {
    m_ready = false;
    if (filenames.size() < 3) {
      m_numGroups = 0;
      m_startSubjects = 0;
      m_startSNPs = 0;
      m_startDosages = 0;
      m_usesTempFiles = false;
      if (m_outfile.good()) {
        if (filenames.size() == 2) {
          if (filenames[1] != "") {
            m_usesTempFiles = true;
            m_chromosomeFilename = filenames[1] + std::string(".chr");
            m_SNPFilename = filenames[1] + std::string(".snp");
            m_refFilename = filenames[1] + std::string(".ref");
            m_altFilename = filenames[1] + std::string(".alt");
            m_bpFilename = filenames[1] + std::string(".bp");
            m_altFreqFilename = filenames[1] + std::string(".altFreq");
            m_mafFilename = filenames[1] + std::string(".maf");
            m_avgCallFilename = filenames[1] + std::string(".avgCall");
            m_rSqFilename = filenames[1] + std::string(".rsq");
            m_chromosomeFile.open(m_chromosomeFilename.c_str());
            m_SNPFile.open(m_SNPFilename.c_str());
            m_refFile.open(m_refFilename.c_str());
            m_altFile.open(m_altFilename.c_str());
            m_bpFile.open(m_bpFilename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
            m_altFreqFile.open(m_altFreqFilename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
            m_mafFile.open(m_mafFilename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
            m_avgCallFile.open(m_avgCallFilename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
            m_rSqFile.open(m_rSqFilename.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
            m_ready = m_chromosomeFile.good() && m_SNPFile.good() && m_refFile.good() && m_altFile.good() && m_bpFile.good();
            m_ready = m_ready && m_altFreqFile.good() && m_mafFile.good() && m_avgCallFile.good() && m_rSqFile.good();
          }
        } else {
          m_ready = true;
        }
      }
    }
  }
}

CWriteBinaryDosage4x::~CWriteBinaryDosage4x() {
  m_chromosomeFile.close();
  m_SNPFile.close();
  m_refFile.close();
  m_altFile.close();
  m_bpFile.close();
  m_altFreqFile.close();
  m_mafFile.close();
  m_avgCallFile.close();
  m_rSqFile.close();
}

int CWriteBinaryDosage4x::AddToStringVector(std::vector<std::string> &addToVector, const std::string &stringToAdd) {
  if (stringToAdd == "") {
    if (addToVector.size() != 0) {
      m_ready = false;
      return 1;
    }
    return 0;
  }
  addToVector.push_back(stringToAdd);
  return 0;
}

int CWriteBinaryDosage4x::AddToDoubleVector(std::vector<std::vector<double> > &addToVector,
                                            const std::vector<double> &vectorToAdd) {
  if (vectorToAdd.size() == 0) {
    if (addToVector.size() != 0) {
      m_ready = false;
      return 1;
    }
    return 0;
  }
  if (vectorToAdd.size() != (unsigned int)m_numGroups)
    return 1;
  addToVector.push_back(vectorToAdd);
  return 0;
}

int CWriteBinaryDosage4x::WriteString(const std::vector<std::string> &stringToWrite) {
  std::vector<std::string>::const_iterator stringIt;
  std::ostringstream oss;
  const char zero = 0x0;

  oss.str("");
  for (stringIt = stringToWrite.begin(); stringIt != stringToWrite.end(); ++stringIt)
    oss << *stringIt << '\t';
  m_outfile.write(oss.str().c_str(), oss.str().length() - 1);
  m_outfile.write((char *)&zero, 1);

  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteHeader() {
  const int zero = 0;

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  for (int i = 0; i < 8; ++i)
    m_outfile.write((char *)&zero, sizeof(int));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteGroups(const std::vector<int> &groupSize) {
  std::streampos startSubPos;

  if (!m_ready)
    return 1;
  m_numGroups = groupSize.size();
  if (m_numGroups == 0) {
    m_ready = false;
    return 1;
  }
  m_outfile.seekp((int)Header4pos::numGroups);
  m_outfile.write((char *)&m_numGroups, sizeof(int));
  m_outfile.seekp((int)Header4pos::startGroups);
  m_outfile.write((char *)groupSize.data(), m_numGroups * sizeof(int));
  startSubPos = m_outfile.tellp();
  m_startSubjects = startSubPos;
  m_outfile.seekp((int)Header4pos::startSub);
  m_outfile.write((char *)&m_startSubjects, sizeof(int));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }

  return 0;
}

int CWriteBinaryDosage4x::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  std::streampos endSub, endFam;
  int numSub;
  int subSize, famSize;
  const int zero = 0;

  if (!m_ready)
    return 1;

  numSub = IID.size();
  if (numSub == 0) {
    m_ready = false;
    return 1;
  }
  m_outfile.seekp((int)Header4pos::numSub);
  m_outfile.write((char *)&numSub, sizeof(int));

  m_outfile.seekp(m_startSubjects);
  m_outfile.write((char *)&zero, sizeof(int));
  m_outfile.write((char *)&zero, sizeof(int));
  WriteString(IID);

  endSub = m_outfile.tellp();
  subSize = (int)endSub - m_startSubjects - 2*sizeof(int);
  m_outfile.seekp(m_startSubjects);
  m_outfile.write((char *)&subSize, sizeof(int));

  m_outfile.seekp(endSub);
  if (FID.size() > 0) {
    WriteString(FID);
    endFam = m_outfile.tellp();
    famSize = (int)endFam - (int)endSub;
    m_outfile.seekp(m_startSubjects + sizeof(int));
    m_outfile.write((char *)&famSize, sizeof(int));
    m_outfile.seekp(endSub);
  }

  m_startSNPs = m_outfile.tellp();
  m_outfile.seekp((int)Header4pos::startSNP);
  m_outfile.write((char *)&m_startSNPs, sizeof(int));
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteStringToTempFile(const std::string &stringToWrite, std::fstream &outfile) {
  if (stringToWrite == "") {
    if (outfile.tellp() != 0) {
      m_ready = false;
      return 1;
    }
  }
  outfile << stringToWrite << '\t';
  if (!outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteIntToTempFile(const int valueToWrite, std::fstream &outfile) {
  outfile.write((char *)&valueToWrite, sizeof(int));
  if (!outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteDoubleVectorToTempFile(const std::vector<double> &vectorToWrite, std::fstream &outfile) {
  if (vectorToWrite.size() == 0) {
    if (outfile.tellp() != 0) {
      m_ready = false;
      return 1;
    }
  }
  outfile.write((char *)vectorToWrite.data(), vectorToWrite.size() * sizeof(double));
  if (!outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}
int CWriteBinaryDosage4x::AddSNP(const std::string &chromosome, const std::string &snpID, int location,
                                   const std::string &refAllele, const std::string &altAllele,
                                   const std::vector<double> &altFreq, const std::vector<double> &maf,
                                   const std::vector<double> &avgCall, const std::vector<double> &rSq) {
  if (!m_ready)
    return 1;
  if (chromosome == "") {
    m_ready = false;
    return 1;
  }
  if (m_usesTempFiles) {
    if (WriteStringToTempFile(chromosome, m_chromosomeFile))
      return 1;
    if (WriteStringToTempFile(snpID, m_SNPFile))
      return 1;
    if (WriteStringToTempFile(refAllele, m_refFile))
      return 1;
    if (WriteStringToTempFile(altAllele, m_altFile))
      return 1;
    if (WriteIntToTempFile(location, m_bpFile))
      return 1;
    if (WriteDoubleVectorToTempFile(altFreq, m_altFreqFile))
      return 1;
    if (WriteDoubleVectorToTempFile(maf, m_mafFile))
      return 1;
    if (WriteDoubleVectorToTempFile(avgCall, m_avgCallFile))
      return 1;
    if (WriteDoubleVectorToTempFile(rSq, m_rSqFile))
      return 1;
  } else {
    m_chromosome.push_back(chromosome);
    m_bp.push_back(location);

    if (AddToStringVector(m_snpID, snpID))
      return 1;
    if (AddToStringVector(m_refAllele, refAllele))
      return 1;
    if (AddToStringVector(m_altAllele, altAllele))
      return 1;
    if (AddToDoubleVector(m_altFreq, altFreq))
      return 1;
    if (AddToDoubleVector(m_maf, maf))
      return 1;
    if (AddToDoubleVector(m_avgCall, avgCall))
      return 1;
    if (AddToDoubleVector(m_rSq, rSq))
      return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::GetSNPOptions() {
  // Chrosome and base pair location are now required
  int snpOptions = 0x0014;
  std::vector<std::string>::const_iterator stringIt;
  std::vector<std::vector<double> >::const_iterator vdIt;
  std::string chromosome;
  unsigned int reqSize;

  if (m_snpID.size() != 0)
    snpOptions |= 0x0002;

  reqSize = m_chromosome.size();
  Rcpp::Rcout << "reqSize\t" << reqSize << std::endl;
  if (reqSize == 0)
    return 0;
  stringIt = m_chromosome.begin();
  chromosome = *stringIt;
  ++stringIt;
  for (;stringIt != m_chromosome.end(); ++stringIt) {
    if (*stringIt != chromosome)
      break;
  }
  if (stringIt == m_chromosome.end())
    snpOptions |= 0x0008;

  Rcpp::Rcout << "bpSize\t" << m_bp.size() << std::endl;
  if (m_bp.size() != reqSize)
    return 0;

  Rcpp::Rcout << "Stop 1" << std::endl;
  if (m_refAllele.size() != 0) {
    if (m_refAllele.size() != reqSize)
      return 0;
    if (m_altAllele.size() != reqSize)
      return 0;
    snpOptions |= 0x0060;
  }

  Rcpp::Rcout << "Stop 2" << std::endl;
  if (m_altFreq.size() != 0) {
    if (m_altFreq.size() != m_numGroups)
      return 0;
    for (vdIt = m_altFreq.begin(); vdIt != m_altFreq.end(); ++vdIt) {
      if (vdIt->size() != reqSize)
        return 0;
    }
    snpOptions |= 0x0080;
  }

  Rcpp::Rcout << "Stop 3" << std::endl;
  if (m_maf.size() != 0) {
    if (m_maf.size() != m_numGroups)
      return 0;
    for (vdIt = m_maf.begin(); vdIt != m_maf.end(); ++vdIt) {
      if (vdIt->size() != reqSize)
        return 0;
    }
    snpOptions |= 0x0100;
  }

  Rcpp::Rcout << "Stop 4" << std::endl;
  if (m_avgCall.size() != 0) {
    if (m_avgCall.size() != m_numGroups)
      return 0;
    for (vdIt = m_avgCall.begin(); vdIt != m_avgCall.end(); ++vdIt) {
      if (vdIt->size() != reqSize)
        return 0;
    }
    snpOptions |= 0x0200;
  }

  Rcpp::Rcout << "Stop 5" << std::endl;
  if (m_rSq.size() != 0) {
    if (m_rSq.size() != m_numGroups)
      return 0;
    for (vdIt = m_rSq.begin(); vdIt != m_rSq.end(); ++vdIt) {
      if (vdIt->size() != reqSize)
        return 0;
    }
    snpOptions |= 0x0400;
  }
  return snpOptions;
}

int CWriteBinaryDosage4x::WriteStringVectorToFile(const std::vector<std::string> &stringToWrite, int sizeLocation) {
  std::streampos startPos, endPos;
  int stringSize;

  startPos = m_outfile.tellp();
  WriteString(stringToWrite);
  endPos = m_outfile.tellp();
  stringSize = (int)(endPos - startPos);
  m_outfile.seekp(sizeLocation);
  m_outfile.write((char *)&stringSize, sizeof(int));
  m_outfile.seekp(endPos);
  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }
  return 0;
}

int CWriteBinaryDosage4x::WriteSNPs() {
  const int zero = 0;
  int numSNPs;
  int snpOptions;
  std::vector<std::string> singleChromosome;
  std::vector<std::vector<double> >::iterator dvIt;
  std::vector<double>::iterator dIt;

  if (!m_ready)
    return 1;

  numSNPs = m_chromosome.size();
  if (numSNPs == 0) {
    m_ready = false;
    return 1;
  }
  m_outfile.seekp((int)Header4pos::numSNPs);
  m_outfile.write((char *)&numSNPs, sizeof(int));
  snpOptions = GetSNPOptions();
  Rcpp::Rcout << "SNP options\t" << std::hex << snpOptions << std::dec << std::endl;
  if (snpOptions == 0) {
    m_ready = false;
    return 1;
  }
  m_outfile.seekp((int)Header4pos::snpOptions);
  m_outfile.write((char *)&snpOptions, sizeof(int));
  m_outfile.seekp(m_startSNPs);
  m_outfile.write((char *)&zero, sizeof(int));
  m_outfile.write((char *)&zero, sizeof(int));
  m_outfile.write((char *)&zero, sizeof(int));
  m_outfile.write((char *)&zero, sizeof(int));
  if (snpOptions & 0x0002)
    WriteStringVectorToFile(m_snpID, m_startSNPs);
  if (snpOptions & 0x0004) {
    if (snpOptions & 0x0008) {
      singleChromosome.push_back(m_chromosome[0]);
      WriteStringVectorToFile(singleChromosome, m_startSNPs + sizeof(int));
    } else {
      WriteStringVectorToFile(m_chromosome, m_startSNPs + sizeof(int));
    }
  }
  m_outfile.write((char *)m_bp.data(), m_bp.size() * sizeof(int));
  if (snpOptions | 0x0020) {
    WriteStringVectorToFile(m_refAllele, m_startSNPs + 2 * sizeof(int));
    WriteStringVectorToFile(m_altAllele, m_startSNPs + 3 * sizeof(int));
  }

  if (snpOptions | 0x0080) {
    for (dvIt = m_altFreq.begin(); dvIt != m_altFreq.end(); ++dvIt) {
      for (dIt = dvIt->begin(); dIt != dvIt->end(); ++dIt) {
        Rcpp::Rcout << "Allele freq\t" << (*dvIt)[0] << std::endl;
        m_outfile.write((char *)dvIt->data(), m_numGroups * sizeof(double));
      }
    }
  }
  if (snpOptions | 0x0100) {
    for (dvIt = m_maf.begin(); dvIt != m_maf.end(); ++dvIt)
      m_outfile.write((char *)dvIt->data(), m_numGroups * sizeof(double));
  }
  if (snpOptions | 0x0200) {
    for (dvIt = m_avgCall.begin(); dvIt != m_avgCall.end(); ++dvIt)
      m_outfile.write((char *)dvIt->data(), m_numGroups * sizeof(double));
  }
  if (snpOptions | 0x0400) {
    for (dvIt = m_rSq.begin(); dvIt != m_rSq.end(); ++dvIt)
      m_outfile.write((char *)dvIt->data(), m_numGroups * sizeof(double));
  }

  m_startDosages = (int)m_outfile.tellp();
  m_outfile.seekp((int)Header4pos::startDosage);
  m_outfile.write((char *)&m_startDosages, sizeof(int));
  m_outfile.seekp(m_startDosages);

  if (!m_outfile.good()) {
    m_ready = false;
    return 1;
  }

  return 0;
}

int CWriteBinaryDosage4x::WriteAllSNPs(const std::vector<std::string> &chromosome,
                                              const std::vector<std::string> &snpID,
                                              const std::vector<int> bp,
                                              const std::vector<std::string> &refAllele,
                                              const std::vector<std::string> &altAllele,
                                              const std::vector<std::vector<double> > &altFreq,
                                              const std::vector<std::vector<double> > &maf,
                                              const std::vector<std::vector<double> > &avgCall,
                                              const std::vector<std::vector<double> > &rSq) {
  int numSNPs;
  std::vector<std::string>::iterator chrIt;
  std::vector<std::vector<double > >::const_iterator vdIt;

  Rcpp::Rcout << chromosome.size() << '\t' << snpID.size() << '\t' << bp.size() << '\t'
              << refAllele.size() << '\t' << altAllele.size() << '\t' << altFreq.size() << '\t'
              << maf.size() << '\t' << avgCall.size() << '\t' << rSq.size() << '\t' << m_numGroups << std::endl;


  if (bp.size() == 0 || chromosome.size() == 0) {
    m_ready = false;
    return 1;
  }
  if (chromosome.size() != 1 && chromosome.size() != bp.size()) {
    m_ready = false;
    return 1;
  }
  if (snpID.size() != 0 && snpID.size() != bp.size()) {
    m_ready = false;
    return 1;
  }
  if (refAllele.size() != 0 && refAllele.size() != bp.size()) {
    m_ready = false;
    return 1;
  }
  if (altAllele.size() != 0 && altAllele.size() != bp.size()) {
    m_ready = false;
    return 1;
  }
  if (altFreq.size() != 0 && altFreq.size() != (unsigned)m_numGroups) {
    m_ready = false;
    return 1;
  }
  for (vdIt = altFreq.begin(); vdIt != altFreq.end(); ++vdIt) {
    if (vdIt->size() != bp.size()) {
      m_ready = false;
      return 1;
    }
  }
  if (maf.size() != 0 && maf.size() != (unsigned)m_numGroups) {
    m_ready = false;
    return 1;
  }
  for (vdIt = maf.begin(); vdIt != maf.end(); ++vdIt) {
    if (vdIt->size() != bp.size()) {
      m_ready = false;
      return 1;
    }
  }
  if (avgCall.size() != 0 && avgCall.size() != (unsigned)m_numGroups) {
    m_ready = false;
    return 1;
  }
  for (vdIt = avgCall.begin(); vdIt != avgCall.end(); ++vdIt) {
    if (vdIt->size() != bp.size()) {
      m_ready = false;
      return 1;
    }
  }
  if (rSq.size() != 0 && rSq.size() != (unsigned)m_numGroups) {
    m_ready = false;
    return 1;
  }
  for (vdIt = rSq.begin(); vdIt != rSq.end(); ++vdIt) {
    if (vdIt->size() != bp.size()) {
      m_ready = false;
      return 1;
    }
  }

  numSNPs = bp.size();
  if (chromosome.size() == 1) {
    m_chromosome.resize(numSNPs);
    std::fill(m_chromosome.begin(), m_chromosome.end(), chromosome[0]);
  } else {
    m_chromosome = chromosome;
  }
  m_snpID = snpID;
  m_bp = bp;
  m_refAllele = refAllele;
  m_altAllele = altAllele;
  m_altFreq = altFreq;
  m_maf = maf;
  m_avgCall = avgCall;
  m_rSq = rSq;

  Rcpp::Rcout << "Before WriteSNPs" << std::endl;
  return WriteSNPs();
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage41
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage41::CWriteBinaryDosage41(const std::vector<std::string> &filenames) : CWriteBinaryDosage4x(filenames) {
  m_version[1] = 0x04;
  m_version[3] = 0x01;
}

int CWriteBinaryDosage41::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddDosagesOnly(geneticValues, 10000.);
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage42
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage42::CWriteBinaryDosage42(const std::vector<std::string> &filenames) : CWriteBinaryDosage4x(filenames) {
  m_version[1] = 0x04;
  m_version[3] = 0x02;
}

int CWriteBinaryDosage42::AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) {
  return AddGeneticValues32(geneticValues);
}
///////////////////////////////////////////////////////////////////////////////
//                            Test code
///////////////////////////////////////////////////////////////////////////////
int CreateSubjectIDs(std::vector<std::string> &subID, const int startID) {
  int i;
  std::string iid;

  for (i = 1; i < 6; ++i) {
    iid = "Subject_" + std::to_string(startID + i);
    subID.push_back(iid);
  }
  return 0;
}

int CreateSNPs(std::vector<std::string> &snpID, std::vector<std::string> &chromosome,
               std::vector<int> &bp, std::vector<std::string> &refAllele, std::vector<std::string> &altAllele) {
  snpID.push_back("SNP1");
  snpID.push_back("SNP2");
  snpID.push_back("SNP3");
  chromosome.push_back("1");
  chromosome.push_back("1");
  chromosome.push_back("1");
  bp.push_back(1001);
  bp.push_back(2001);
  bp.push_back(3001);
  refAllele.push_back("C");
  refAllele.push_back("A");
  refAllele.push_back("C");
  altAllele.push_back("T");
  altAllele.push_back("T");
  altAllele.push_back("G");
  return 0;
}

int CreateExtraValues(std::vector<std::vector<double> > &altFreq, std::vector<std::vector<double> > &maf,
                      std::vector<std::vector<double> > &avgCall, std::vector<std::vector<double> > &rSq) {
  altFreq.resize(3);
  maf.resize(3);
  avgCall.resize(3);
  rSq.resize(3);
  altFreq[0].push_back(0.8);
  altFreq[1].push_back(0.2);
  altFreq[2].push_back(0.1);
  maf[0].push_back(0.2);
  maf[1].push_back(0.2);
  maf[2].push_back(0.1);
  return 0;
}

int CreateGeneticValues(std::vector<std::vector<std::vector<double> > > &geneticValues) {
  const double gData[60] = {
    0., 1., 0., 0.,
    1., 0., 1., 0.,
    2., 0., 0., 1.,
    0.5, 0.5, 0.5, 0.,
    1.5, 0., 0.5, 0.5,
    0.1, 0.9, 0.1, 0.,
    0.2, 0.79, 0.2, 0.,
    0.3, 0.71, 0.28, 0.01,
    0., 1., 0., 0.,
    1.2, 0., 0.8, 0.2,
    1.5, 0., 0.5, 0.5,
    0.5, 0.5, 0.5, 0.,
    2., 0., 0., 1.,
    1., 0., 1., 0.,
    0., 1., 0., 0.
  };
  std::vector<std::vector<std::vector<double > > >::iterator v3It;
  std::vector<std::vector<double > >::iterator v2It;
  std::vector<double>::iterator dIt;
  const double *d1, *d2, *d3;

  geneticValues.resize(3);
  d1 = gData;
  for (v3It = geneticValues.begin(); v3It != geneticValues.end(); ++v3It, d1 += 20) {
    d2 = d1;
    v3It->resize(4);
    for (v2It = v3It->begin(); v2It != v3It->end(); ++v2It, ++d2) {
      d3 = d2;
      v2It->resize(5);
      for (dIt = v2It->begin(); dIt != v2It->end(); ++dIt, d3 += 4) {
        *dIt = *d3;
      }
    }
  }
  return 0;
}

int TestWriteBD(CWriteBinaryDosage *bdf, const int startID) {
  std::vector<int> groups;
  std::vector<std::string> fid, iid;
  std::vector<std::string> snpID, chromosome, refAllele, altAllele;
  std::vector<int> bp;
  std::vector<std::vector<double> > altFreq, maf, avgCall, rSq;
  std::vector<std::vector<std::vector<double> > > geneticValues;

  bdf->WriteHeader();
  groups.push_back(5);
  bdf->WriteGroups(groups);
  CreateSubjectIDs(iid, startID);
  bdf->WriteSubjects(fid, iid);
  CreateSNPs(snpID, chromosome, bp, refAllele, altAllele);
  CreateExtraValues(altFreq, maf, avgCall, rSq);
  for (int i = 0; i < 3; ++i)
    bdf->AddSNP(chromosome[i], snpID[i], bp[i], refAllele[i], altAllele[i], altFreq[i], maf[i], avgCall[i], rSq[i]);
  bdf->WriteSNPs();
  CreateGeneticValues(geneticValues);
  for (int i = 0; i < 3; ++i)
    bdf->AddGeneticValues(geneticValues[i]);
  return 0;
}
//' Function to test writing of binary dosage files
//'
//' Function to test writing of binary dosage files
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
int TestWriteBinaryDosage() {
  std::vector<std::string> f11Files, f12Files, f21Files, f22Files, f31Files, f32Files, f41Files, f42Files, f42aFiles, f42bFiles, f42cFiles;

  f11Files.push_back("Test/Test2.Format11.bdose");
  f11Files.push_back("Test/Test2.Format11.fam");
  f11Files.push_back("Test/Test2.Format11.map");
  f12Files.push_back("Test/Test2.Format12.bdose");
  f12Files.push_back("Test/Test2.Format12.fam");
  f12Files.push_back("Test/Test2.Format12.map");
  f21Files.push_back("Test/Test2.Format21.bdose");
  f21Files.push_back("Test/Test2.Format21.fam");
  f21Files.push_back("Test/Test2.Format21.map");
  f22Files.push_back("Test/Test2.Format22.bdose");
  f22Files.push_back("Test/Test2.Format22.fam");
  f22Files.push_back("Test/Test2.Format22.map");
  f31Files.push_back("Test/Test2.Format31.bdose");
  f31Files.push_back("Test/Test2.Format31.fam");
  f31Files.push_back("Test/Test2.Format31.map");
  f32Files.push_back("Test/Test2.Format32.bdose");
  f32Files.push_back("Test/Test2.Format32.fam");
  f32Files.push_back("Test/Test2.Format32.map");
  f41Files.push_back("Test/Test2.Format41.bdose");
  f42Files.push_back("Test/Test2.Format42.bdose");
  f42aFiles.push_back("Test/Test2.Format42a.bdose");
  f42bFiles.push_back("Test/Test2.Format42b.bdose");
  f42cFiles.push_back("Test/Test2.Format42c.bdose");

  CWriteBinaryDosage11 f11_1(f11Files);
  CWriteBinaryDosage12 f12_1(f12Files);
  CWriteBinaryDosage21 f21_1(f21Files);
  CWriteBinaryDosage22 f22_1(f22Files);
  CWriteBinaryDosage31 f31_1(f31Files);
  CWriteBinaryDosage32 f32_1(f32Files);
  CWriteBinaryDosage41 f41_1(f41Files);
  CWriteBinaryDosage42 f42_1(f42Files);
  CWriteBinaryDosage42 f42a_1(f42aFiles);
  CWriteBinaryDosage42 f42b_1(f42bFiles);
  CWriteBinaryDosage42 f42c_1(f42cFiles);

  TestWriteBD(&f11_1, 1);
  TestWriteBD(&f12_1, 1);
  TestWriteBD(&f21_1, 1);
  TestWriteBD(&f22_1, 1);
  TestWriteBD(&f31_1, 1);
  TestWriteBD(&f32_1, 1);
  TestWriteBD(&f41_1, 1);
  TestWriteBD(&f42_1, 1);
  TestWriteBD(&f42a_1, 6);
  TestWriteBD(&f42b_1, 11);
  TestWriteBD(&f42c_1, 16);

  return 0;
}
