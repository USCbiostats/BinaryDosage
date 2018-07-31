#include <iostream>
#include <fstream>
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

int CWriteMultifileBinaryDosage::WriteSNP(const std::string &chromosome, const std::string &snpID, int location,
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
  int numGroups;
  std::streampos startSubPos;
  int startSubjects;

  numGroups = groupSize.size();
  m_outfile.seekp((int)Header4pos::numGroups);
  m_outfile.write((char *)&numGroups, sizeof(int));
  m_outfile.seekp((int)Header4pos::startGroups);
  m_outfile.write((char *)groupSize.data(), numGroups * sizeof(int));
  startSubPos = m_outfile.tellp();
  startSubjects = startSubPos;
  m_outfile.seekp((int)Header4pos::startSub);
  m_outfile.write((char *)&startSubjects, sizeof(int));
  if (!m_outfile.good())
    return 1;

  return 0;
}

int CWriteBinaryDosage4x::WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) {

  return 0;
}

int CWriteBinaryDosage4x::WriteSNP(const std::string &chromosome, const std::string &snpID, int location,
                                   const std::string &refAllele, const std::string &altAllele) {
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

int TestWriteBD(CWriteBinaryDosage *bdf) {
  std::vector<int> groups;
  std::vector<std::string> fid, iid;
  std::vector<std::string> snpID, chromosome, refAllele, altAllele;
  std::vector<int> bp;

  bdf->WriteHeader();
  groups.push_back(5);
  bdf->WriteGroups(groups);
  CreateSubjectIDs(iid);
  bdf->WriteSubjects(fid, iid);
  CreateSNPs(snpID, chromosome, bp, refAllele, altAllele);
  for (int i = 0; i < 3; ++i)
    bdf->WriteSNP(chromosome[i], snpID[i], bp[i], refAllele[i], altAllele[i]);
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
