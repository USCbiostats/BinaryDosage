#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "ReadBinaryDosage.h"
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
//                            CReadBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosageX::CReadBinaryDosageX(const std::string &filename) {
  m_filename = filename;
  memset(m_version, 0, sizeof(m_version));
  m_mainVersion = 0;
  m_subVersion = 0;
  m_numSubjects = 0;
  m_numSNPs = 0;
  m_numGroups = 0;
  m_subjectOptions = 0;
  m_snpOptions = 0;
  m_usesFamilyID = false;
}

CReadBinaryDosageX::~CReadBinaryDosageX() {
  m_infile.close();
}

int CReadBinaryDosageX::ReadVersion() {
  m_infile.read((char *)m_version, 4);
  if (!m_infile.good())
    return 1;
  return 0;
}

int CReadBinaryDosageX::ReadHeader() {
  const char bdHeader[4] = { 'b', 'o', 's', 'e' };
  char header[4];
  if (!m_infile.is_open())
    return 1;

  m_infile.read(header, 4);
  if (!m_infile.good())
    return 1;
  if (memcmp(bdHeader, header, 4))
    return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadMultifileBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CReadMultifileBinaryDosage::CReadMultifileBinaryDosage(const std::string &filename) : CReadBinaryDosageX(filename) {
  std::string outputFilename;

  outputFilename = m_filename + std::string(".fam");
  m_famFile.open(outputFilename.c_str());
  outputFilename = m_filename + std::string(".map");
  m_mapFile.open(outputFilename.c_str());
  outputFilename = m_filename + std::string(".bdose");
  m_infile.open(outputFilename.c_str(), std::ios_base::in | std::ios_base::binary);
}

CReadMultifileBinaryDosage::~CReadMultifileBinaryDosage() {
  m_famFile.close();
  m_mapFile.close();
}

int CReadMultifileBinaryDosage::ReadSubjects() {
  std::vector<std::string>::const_iterator iid, fid;
  std::string junk, fam, sub;
  std::istringstream iss;
  int numCol;


  if (!m_famFile.is_open())
    return 1;

  std::getline(m_famFile, junk);
  if (m_famFile.fail())
    return 1;
  numCol = 0;
  iss.str(junk);

  if (FID.size() == 0) {
    for (std::vector<std::string>::const_iterator iid = IID.begin(); iid != IID.end(); ++iid)
      m_famFile << *iid << "\t0\t0\t9\t9" << std::endl;
  } else {
    for (iid = IID.begin(), fid = FID.begin(); iid != IID.end(); ++iid, ++fid)
      m_famFile << *fid << '\t' << *iid << "\t0\t0\t9\t9" << std::endl;
  }
  return 0;
}

int CReadMultifileBinaryDosage::ReadSNP(const std::string &chromosome, const std::string &snpID, int location,
                                          const std::string &refAllele, const std::string &altAllele) {
  if (!m_mapFile.good())
    return 1;
  m_mapFile << chromosome << '\t' << snpID << "\t0\t" << location << '\t' << refAllele << '\t' << altAllele << std::endl;
  if (!m_mapFile.good())
    return 1;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage11
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage11::CReadBinaryDosage11(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage11::ReadHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x01};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage12
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage12::CReadBinaryDosage12(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage12::ReadHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x02};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage21
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage21::CReadBinaryDosage21(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage21::ReadHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x01};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage22
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage22::CReadBinaryDosage22(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage22::ReadHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x02};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage31
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage31::CReadBinaryDosage31(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage31::ReadHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x01};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

int CReadBinaryDosage31::ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  numSubjects = IID.size();
  m_infile.Read((char *)&numSubjects, sizeof(int));
  if (!m_infile.good())
    return 1;
  return CReadMultifileBinaryDosage::ReadSubjects(FID, IID);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage32
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage32::CReadBinaryDosage32(const std::string &filename) : CReadMultifileBinaryDosage(filename) {}

int CReadBinaryDosage32::ReadHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x02};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  return ReadVersion(version);
}

int CReadBinaryDosage32::ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {
  int numSubjects;

  numSubjects = IID.size();
  m_infile.Read((char *)&numSubjects, sizeof(int));
  if (!m_infile.good())
    return 1;
  return CReadMultifileBinaryDosage::ReadSubjects(FID, IID);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage4x
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage4x::CReadBinaryDosage4x(const std::string &filename) : CReadBinaryDosage(filename) {
  std::string infilename;

  infilename = m_filename + std::string(".bdose");
  m_infile.open(infilename.c_str(), std::ios_base::out | std::ios_base::binary);
}

int CReadBinaryDosage4x::ReadHeader() {
  const int zero = 0;

  for (int i = 0; i < 8; ++i)
    m_infile.Read((char *)&zero, sizeof(int));
  if (!m_infile.good())
    return 1;
  return 0;
}

int CReadBinaryDosage4x::ReadGroups(const std::vector<int> &groupSize) {
  int numGroups;
  std::streampos startSubPos;
  int startSubjects;

  numGroups = groupSize.size();
  m_infile.seekp((int)Header4pos::numGroups);
  m_infile.Read((char *)&numGroups, sizeof(int));
  m_infile.seekp((int)Header4pos::startGroups);
  m_infile.Read((char *)groupSize.data(), numGroups * sizeof(int));
  startSubPos = m_infile.tellp();
  startSubjects = startSubPos;
  m_infile.seekp((int)Header4pos::startSub);
  m_infile.Read((char *)&startSubjects, sizeof(int));
  if (!m_infile.good())
    return 1;

  return 0;
}

int CReadBinaryDosage4x::ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {

  return 0;
}

int CReadBinaryDosage4x::ReadSNP(const std::string &chromosome, const std::string &snpID, int location,
                                   const std::string &refAllele, const std::string &altAllele) {
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage41
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage41::CReadBinaryDosage41(const std::string &filename) : CReadBinaryDosage4x(filename) {}

int CReadBinaryDosage41::ReadHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x01};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  if (ReadVersion(version))
    return 1;
  return CReadBinaryDosage4x::ReadHeader();
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage42
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage42::CReadBinaryDosage42(const std::string &filename) : CReadBinaryDosage4x(filename) {}

int CReadBinaryDosage42::ReadHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x02};

  if (CReadBinaryDosage::ReadHeader())
    return 1;
  if (ReadVersion(version))
    return 1;
  return CReadBinaryDosage4x::ReadHeader();
}

///////////////////////////////////////////////////////////////////////////////
//                            Test code
///////////////////////////////////////////////////////////////////////////////

int TestReadBD(CReadBinaryDosage *bdf) {
  std::vector<int> groups;
  std::vector<std::string> fid, iid;
  std::vector<std::string> snpID, chromosome, refAllele, altAllele;
  std::vector<int> bp;

  bdf->ReadHeader();
  groups.push_back(5);
  bdf->ReadGroups(groups);
  CreateSubjectIDs(iid);
  bdf->ReadSubjects(fid, iid);
  CreateSNPs(snpID, chromosome, bp, refAllele, altAllele);
  for (int i = 0; i < 3; ++i)
    bdf->ReadSNP(chromosome[i], snpID[i], bp[i], refAllele[i], altAllele[i]);
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
int TestReadBinaryDosage() {
  CReadBinaryDosage11 f11("Test/Test.Format11");
  CReadBinaryDosage12 f12("Test/Test.Format12");
  CReadBinaryDosage21 f21("Test/Test.Format21");
  CReadBinaryDosage22 f22("Test/Test.Format22");
  CReadBinaryDosage31 f31("Test/Test.Format31");
  CReadBinaryDosage32 f32("Test/Test.Format32");
  CReadBinaryDosage41 f41("Test/Test.Format41");
  CReadBinaryDosage42 f42("Test/Test.Format42");

  TestReadBD(&f11);
  TestReadBD(&f12);
  TestReadBD(&f21);
  TestReadBD(&f22);
  TestReadBD(&f31);
  TestReadBD(&f32);
  TestReadBD(&f41);
  TestReadBD(&f42);

  return 0;
}
