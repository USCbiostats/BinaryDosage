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

CWriteBinaryDosage::CWriteBinaryDosage(const std::string &filename) {
  m_filename = filename;
}

CWriteBinaryDosage::~CWriteBinaryDosage() {
  m_outfile.close();
}

int CWriteBinaryDosage::WriteVersion(const char *version) {
  m_outfile.write(version, 4);
  if (!m_outfile.good())
    return 1;
  return 0;
}

int CWriteBinaryDosage::WriteHeader() {
  const char header[4] = { 'b', 'o', 's', 'e' };
  if (!m_outfile.is_open())
    return 1;

  m_outfile.write(header, 4);
  if (!m_outfile.good())
    return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteMultifileBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CWriteMultifileBinaryDosage::CWriteMultifileBinaryDosage(const std::string &filename) : CWriteBinaryDosage(filename) {
  std::string outputFilename;

  outputFilename = m_filename + std::string(".fam");
  m_famFile.open(outputFilename.c_str());
  outputFilename = m_filename + std::string(".map");
  m_mapFile.open(outputFilename.c_str());
  outputFilename = m_filename + std::string(".bdose");
  m_outfile.open(outputFilename.c_str(), std::ios_base::out | std::ios_base::binary);
}

CWriteMultifileBinaryDosage::~CWriteMultifileBinaryDosage() {
  m_famFile.close();
  m_mapFile.close();
}

int CWriteMultifileBinaryDosage::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  std::vector<std::string>::const_iterator iid, fid;

  if (!m_famFile.is_open())
    return 1;

  if (FID.size() == 0) {
    for (std::vector<std::string>::const_iterator iid = IID.begin(); iid != IID.end(); ++iid)
      m_famFile << *iid << "\t0\t0\t9\t9" << std::endl;
  } else {
    for (iid = IID.begin(), fid = FID.begin(); iid != IID.end(); ++iid, ++fid)
      m_famFile << *fid << '\t' << *iid << "\t0\t0\t9\t9" << std::endl;
  }
  return 0;
}

int CWriteMultifileBinaryDosage::AddSNP(const std::string &chromosome, const std::string &snpID, int location,
                                          const std::string &refAllele, const std::string &altAllele,
                                          const std::vector<double> &altFreq, const std::vector<double> &maf,
                                          const std::vector<double> &avgCall, const std::vector<double> &rSq) {
  if (!m_mapFile.good())
    return 1;
  m_mapFile << chromosome << '\t' << snpID << "\t0\t" << location << '\t' << refAllele << '\t' << altAllele << std::endl;
  if (!m_mapFile.good())
    return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage11
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage11::CWriteBinaryDosage11(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage11::WriteHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x01};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage12
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage12::CWriteBinaryDosage12(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage12::WriteHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x02};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage21
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage21::CWriteBinaryDosage21(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage21::WriteHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x01};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage22
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage22::CWriteBinaryDosage22(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage22::WriteHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x02};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage31
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage31::CWriteBinaryDosage31(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage31::WriteHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x01};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

int CWriteBinaryDosage31::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  numSubjects = IID.size();
  m_outfile.write((char *)&numSubjects, sizeof(int));
  if (!m_outfile.good())
    return 1;
  return CWriteMultifileBinaryDosage::WriteSubjects(FID, IID);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage32
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage32::CWriteBinaryDosage32(const std::string &filename) : CWriteMultifileBinaryDosage(filename) {}

int CWriteBinaryDosage32::WriteHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x02};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  return WriteVersion(version);
}

int CWriteBinaryDosage32::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  numSubjects = IID.size();
  m_outfile.write((char *)&numSubjects, sizeof(int));
  if (!m_outfile.good())
    return 1;
  return CWriteMultifileBinaryDosage::WriteSubjects(FID, IID);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage4x
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage4x::CWriteBinaryDosage4x(const std::string &filename) : CWriteBinaryDosage(filename) {
  std::string outFilename;

  outFilename = m_filename + std::string(".bdose");
  m_outfile.open(outFilename.c_str(), std::ios_base::out | std::ios_base::binary);
  m_numGroups = 0;
  m_startSubjects = 0;
  m_startSNPs = 0;
  m_startDosages = 0;
}

int CWriteBinaryDosage4x::AddToStringVector(std::vector<std::string> &addToVector, const std::string &stringToAdd) {
  if (stringToAdd == "") {
    if (addToVector.size() != 0)
      return 1;
    return 0;
  }
  addToVector.push_back(stringToAdd);
  return 0;
}

int CWriteBinaryDosage4x::AddToDoubleVector(std::vector<std::vector<double> > &addToVector,
                                            const std::vector<double> &vectorToAdd) {
  if (vectorToAdd.size() == 0) {
    if (addToVector.size() != 0)
      return 1;
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

  if (!m_outfile.good())
    return 1;
  return 0;
}

int CWriteBinaryDosage4x::WriteHeader() {
  const int zero = 0;

  for (int i = 0; i < 8; ++i)
    m_outfile.write((char *)&zero, sizeof(int));
  if (!m_outfile.good())
    return 1;
  return 0;
}

int CWriteBinaryDosage4x::WriteGroups(const std::vector<int> &groupSize) {
  std::streampos startSubPos;

  m_numGroups = groupSize.size();
  m_outfile.seekp((int)Header4pos::numGroups);
  m_outfile.write((char *)&m_numGroups, sizeof(int));
  m_outfile.seekp((int)Header4pos::startGroups);
  m_outfile.write((char *)groupSize.data(), m_numGroups * sizeof(int));
  startSubPos = m_outfile.tellp();
  m_startSubjects = startSubPos;
  m_outfile.seekp((int)Header4pos::startSub);
  m_outfile.write((char *)&m_startSubjects, sizeof(int));
  if (!m_outfile.good())
    return 1;

  return 0;
}

int CWriteBinaryDosage4x::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  std::streampos endSub, endFam;
  int numSub;
  int subSize, famSize;
  const int zero = 0;

  numSub = IID.size();
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
  return 0;
}

int CWriteBinaryDosage4x::AddSNP(const std::string &chromosome, const std::string &snpID, int location,
                                   const std::string &refAllele, const std::string &altAllele,
                                   const std::vector<double> &altFreq, const std::vector<double> &maf,
                                   const std::vector<double> &avgCall, const std::vector<double> &rSq) {
  if (chromosome == "")
    return 1;
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
  return 0;
}

int CWriteBinaryDosage4x::GetSNPOptions() {
  // Chrosome and base pair location are now required
  int snpOptions = 0x0014;
  std::vector<std::string>::const_iterator stringIt;
  std::string chromosome;
  unsigned int reqSize;

  if (m_snpID.size() != 0)
    snpOptions |= 0x0002;

  reqSize = m_chromosome.size();
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

  if (m_bp.size() != reqSize)
    return 0;

  if (m_refAllele.size() != 0) {
    if (m_refAllele.size() != reqSize)
      return 0;
    if (m_altAllele.size() != reqSize)
      return 0;
    snpOptions |= 0x0060;
  }

  if (m_altFreq.size() != 0) {
    if (m_altFreq.size() != reqSize)
      return 0;
    snpOptions |= 0x0080;
  }

  if (m_maf.size() != 0) {
    if (m_maf.size() != reqSize)
      return 0;
    snpOptions |= 0x0100;
  }

  if (m_avgCall.size() != 0) {
    if (m_avgCall.size() != reqSize)
      return 0;
    snpOptions |= 0x0200;
  }

  if (m_rSq.size() != 0) {
    if (m_rSq.size() != reqSize)
      return 0;
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
  return 0;
}

int CWriteBinaryDosage4x::WriteSNPs() {
  const int zero = 0;
  int numSNPs;
  int snpOptions;
  std::vector<std::string> singleChromosome;
  std::vector<std::vector<double> >::iterator dvIt;

  numSNPs = m_chromosome.size();
  if (numSNPs == 0)
    return 1;
  m_outfile.seekp((int)Header4pos::numSNPs);
  m_outfile.write((char *)&numSNPs, sizeof(int));
  m_outfile.seekp((int)Header4pos::startSNP);
  m_outfile.write((char *)&m_startSNPs, sizeof(int));
  snpOptions = GetSNPOptions();
  Rcpp::Rcout << std::hex << snpOptions << std::dec << std::endl;
  if (snpOptions == 0)
    return 1;
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
  m_outfile.write((char *)m_bp.data(), m_bp.size() * sizeof(double));
  if (snpOptions | 0x0020) {
    WriteStringVectorToFile(m_refAllele, m_startSNPs + 2 * sizeof(int));
    WriteStringVectorToFile(m_altAllele, m_startSNPs + 3 * sizeof(int));
  }

  if (snpOptions | 0x0080) {
    for (dvIt = m_altFreq.begin(); dvIt != m_altFreq.end(); ++dvIt)
      m_outfile.write((char *)dvIt->data(), m_numGroups * sizeof(double));
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

  if (m_outfile.fail())
    return 1;

  m_startDosages = (int)m_outfile.tellp();
  m_outfile.seekp((int)Header4pos::startDosage);
  m_outfile.write((char *)&m_startDosages, sizeof(int));
  m_outfile.seekp(m_startDosages);
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage41
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage41::CWriteBinaryDosage41(const std::string &filename) : CWriteBinaryDosage4x(filename) {}

int CWriteBinaryDosage41::WriteHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x01};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  if (WriteVersion(version))
    return 1;
  return CWriteBinaryDosage4x::WriteHeader();
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CWriteBinaryDosage42
//
///////////////////////////////////////////////////////////////////////////////

CWriteBinaryDosage42::CWriteBinaryDosage42(const std::string &filename) : CWriteBinaryDosage4x(filename) {}

int CWriteBinaryDosage42::WriteHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x02};

  if (CWriteBinaryDosage::WriteHeader())
    return 1;
  if (WriteVersion(version))
    return 1;
  return CWriteBinaryDosage4x::WriteHeader();
}

///////////////////////////////////////////////////////////////////////////////
//                            Test code
///////////////////////////////////////////////////////////////////////////////
int CreateSubjectIDs(std::vector<std::string> &subID) {
  int i;
  std::string iid;

  for (i = 1; i < 6; ++i) {
    iid = "Subject_" + std::to_string(i);
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

int TestWriteBD(CWriteBinaryDosage *bdf) {
  std::vector<int> groups;
  std::vector<std::string> fid, iid;
  std::vector<std::string> snpID, chromosome, refAllele, altAllele;
  std::vector<int> bp;
  std::vector<std::vector<double> > altFreq, maf, avgCall, rSq;

  bdf->WriteHeader();
  groups.push_back(5);
  bdf->WriteGroups(groups);
  CreateSubjectIDs(iid);
  bdf->WriteSubjects(fid, iid);
  CreateSNPs(snpID, chromosome, bp, refAllele, altAllele);
  CreateExtraValues(altFreq, maf, avgCall, rSq);
  for (int i = 0; i < 3; ++i)
    bdf->AddSNP(chromosome[i], snpID[i], bp[i], refAllele[i], altAllele[i], altFreq[i], maf[i], avgCall[i], rSq[i]);
  bdf->WriteSNPs();
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
  CWriteBinaryDosage11 f11("Test/Test.Format11");
  CWriteBinaryDosage12 f12("Test/Test.Format12");
  CWriteBinaryDosage21 f21("Test/Test.Format21");
  CWriteBinaryDosage22 f22("Test/Test.Format22");
  CWriteBinaryDosage31 f31("Test/Test.Format31");
  CWriteBinaryDosage32 f32("Test/Test.Format32");
  CWriteBinaryDosage41 f41("Test/Test.Format41");
  CWriteBinaryDosage42 f42("Test/Test.Format42");

  TestWriteBD(&f11);
  TestWriteBD(&f12);
  TestWriteBD(&f21);
  TestWriteBD(&f22);
  TestWriteBD(&f31);
  TestWriteBD(&f32);
  TestWriteBD(&f41);
  TestWriteBD(&f42);

  return 0;
}
