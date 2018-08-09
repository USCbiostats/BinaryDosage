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

CReadBinaryDosageX::CReadBinaryDosageX(const std::vector<std::string> &filenames) {
  m_good = false;
  memset(m_version, 0, sizeof(m_version));
  m_mainVersion = 0;
  m_subVersion = 0;
  m_numSubjects = 0;
  m_numSNPs = 0;
  m_numGroups = 0;
  m_startSubjects = 0;
  m_startSNPs = 0;
  m_startDosages = 8;
  m_subjectOptions = 0;
  m_snpOptions = 0;
  m_currentSNP = -1;
  m_usesFamilyID = false;
  m_scale = 10000.;
  if (filenames.size() == 0) {
    m_filename = "";
  } else {
    m_filename = filenames[0];
    m_infile.open(m_filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (m_infile.good())
      m_good = true;
  }
}

CReadBinaryDosageX::~CReadBinaryDosageX() {
  m_infile.close();
}

int CReadBinaryDosageX::ReadVersion(const char *version) {
  m_infile.read((char *)m_version, 4);
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }
  if (std::memcmp(m_version, version, 4))
    return 1;
  return 0;
}

int CReadBinaryDosageX::ReadHeader() {
  const char bdHeader[4] = { 'b', 'o', 's', 'e' };
  char header[4];

  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }

  m_infile.read(header, 4);
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }
  if (std::memcmp(bdHeader, header, 4)) {
    m_good = false;
    return 1;
  }
  return 0;
}

int CReadBinaryDosageX::ReadGroups() {
  m_numGroups = 1;
  m_groupSize.resize(1);
  m_groupSize[0] = m_numSubjects;
  return 0;
}

int CReadBinaryDosageX::ReadDosageSNP() {
  std::vector<unsigned short>::iterator sIt;
  std::vector<double>::iterator dIt;

  if (!m_good)
    return 1;
  if (m_sReadValues.size() == 0) {
    m_sReadValues.resize(m_numSubjects);
  } else {
    if (m_sReadValues.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_dosage.size() == 0) {
    m_dosage.resize(m_numSubjects);
  } else {
    if (m_dosage.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }

  m_infile.read((char *)m_sReadValues.data(), m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }
  dIt = m_dosage.begin();
  for (sIt = m_sReadValues.begin(); sIt != m_sReadValues.end(); ++sIt, ++dIt)
    *dIt = (double)(*sIt) / m_scale;

  return 0;
}

int CReadBinaryDosageX::ReadFullGeneticSNP() {
  std::vector<unsigned short>::iterator sP1It, sP2It;
  std::vector<double>::iterator dIt, p0It, p1It, p2It;

  if (!m_good)
    return 1;
  if (m_sReadValues.size() == 0) {
    m_sReadValues.resize(2 * m_numSubjects);
  } else {
    if (m_sReadValues.size() != 2 * m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_dosage.size() == 0) {
    m_dosage.resize(m_numSubjects);
  } else {
    if (m_dosage.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p0.size() == 0) {
    m_p0.resize(m_numSubjects);
  } else {
    if (m_p0.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p1.size() == 0) {
    m_p1.resize(m_numSubjects);
  } else {
    if (m_p1.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p2.size() == 0) {
    m_p2.resize(m_numSubjects);
  } else {
    if (m_p2.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }

  m_infile.read((char *)m_sReadValues.data(), 2 * m_numSubjects * sizeof(unsigned short));
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }

  p0It = m_p0.begin();
  p1It = m_p1.begin();
  p2It = m_p2.begin();
  sP1It = m_sReadValues.begin();
  sP2It = sP1It + m_numSubjects;
  for (dIt = m_dosage.begin(); dIt != m_dosage.end(); ++dIt, ++p0It, ++p1It, ++p2It, ++sP1It, ++sP2It) {
    *p1It = (double)(*sP1It) / m_scale;
    *p2It = (double)(*sP2It) / m_scale;
    *p0It = 1. - *p1It - *p2It;
    if (*p0It < 0.)
      *p0It = 0.;
    *dIt = *p1It + *p2It + *p2It;
    if (*dIt > 2.)
      *dIt = 2;
  }

  return 0;
}

int CReadBinaryDosageX::ReadV4GeneticSNP() {
  std::vector<unsigned short>::iterator sP1It, sP2It;
  std::vector<double>::iterator dIt, p0It, p1It, p2It;
  int snpSize;

  if (!m_good)
    return 1;
  if (m_sReadValues.size() == 0) {
    m_sReadValues.resize(4 * m_numSubjects);
  } else {
    if (m_sReadValues.size() != 4 * m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_dosage.size() == 0) {
    m_dosage.resize(m_numSubjects);
  } else {
    if (m_dosage.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p0.size() == 0) {
    m_p0.resize(m_numSubjects);
  } else {
    if (m_p0.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p1.size() == 0) {
    m_p1.resize(m_numSubjects);
  } else {
    if (m_p1.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }
  if (m_p2.size() == 0) {
    m_p2.resize(m_numSubjects);
  } else {
    if (m_p2.size() != m_numSubjects) {
      m_good = false;
      return 1;
    }
  }

  m_infile.read((char *)&snpSize, sizeof(int));
  m_infile.read((char *)m_sReadValues.data(), snpSize);
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }

  p0It = m_p0.begin();
  p1It = m_p1.begin();
  p2It = m_p2.begin();
  sP1It = m_sReadValues.begin();
  sP2It = sP1It + m_numSubjects;
  for (dIt = m_dosage.begin(); dIt != m_dosage.end(); ++dIt, ++p0It, ++p1It, ++p2It, ++sP1It) {
    if ((*sP1It & 0x8000) != 0) {
      *sP1It &= 0x7fff;
      *dIt = (double)(*sP1It) / m_scale;
      if ((*sP2It & 0x7fff) != 0) {
        *sP2It &= 0x7fff;
        *p1It = (double)(*sP2It) / m_scale;
        ++sP2It;
        *p0It = (double)(*sP2It) / m_scale;
        ++sP2It;
        *p2It = (double)(*sP2It) / m_scale;
        ++sP2It;
      } else {
        *p1It = (double)(*sP2It) / m_scale;
        ++sP2It;
        *p2It = (*dIt - *p1It) / 2.;
        *p0It = 1. - *p1It - *p2It;
      }
    } else {
      *dIt = (double)(*sP1It) / m_scale;
      if (*dIt > 1) {
        *p0It = 0.;
        *p2It = *dIt - 1.;
        *p1It = 1. - *p2It;
      } else {
        *p2It = 0.;
        *p1It = *dIt;
        *p0It = 1. - *p1It;
      }
    }
  }

  return 0;
}

void CReadBinaryDosageX::WriteData(std::ostream &outfile) {
  std::vector<double>::iterator dit;

  outfile << "---------------------------------------------------------" << std::endl;
  outfile << "Filename\t\t:\t" << m_filename << std::endl;
  outfile << "Version\t\t\t:\t" << m_mainVersion << '.' << m_subVersion << std::endl;
  outfile << "Number of subjects\t:\t" << m_numSubjects << std::endl;
  outfile << "Subject options\t\t:\t" << std::hex << m_subjectOptions << std::dec << std::endl;
  outfile << "Start subjects\t\t:\t" << m_startSubjects << std::endl;
  outfile << "First subject\t\t:\t";
  if (m_FID.size() > 0)
    outfile << m_FID[0] << '\t';
  if (m_IID.size() > 0)
    outfile << m_IID[0];
  outfile << std::endl;
  outfile << "Last subject\t\t:\t";
  if (m_FID.size() > 0)
    outfile << m_FID[m_FID.size() - 1] << '\t';
  if (m_IID.size() > 0)
    outfile << m_IID[m_IID.size() - 1];
  outfile << std::endl;
  outfile << "Number of groups\t:\t" << m_numGroups << std::endl;
  outfile << "Group sizes\t\t:";
  for (int i = 0; i < m_numGroups; ++i)
    outfile << '\t' << m_groupSize[i];
  outfile << std::endl;
  outfile << "Number of SNPs\t\t:\t" << m_numSNPs << std::endl;
  outfile << "Start of SNPs\t\t:\t" << m_startSNPs << std::endl;
  outfile << "SNP options\t\t:\t" << std::hex << m_snpOptions << std::dec << std::endl;
  outfile << "First SNP\t\t:";
  if (m_SNPID.size() > 0)
    outfile << '\t' << m_SNPID[0];
  if (m_chromosome.size() > 0)
    outfile << '\t' << m_chromosome[0];
  if (m_bp.size() > 0)
    outfile << '\t' << m_bp[0];
  if (m_refAllele.size() > 0)
    outfile << '\t' << m_refAllele[0];
  if (m_altAllele.size() > 0)
    outfile << '\t' << m_altAllele[0];
  outfile << std::endl;
  outfile << "Last SNP\t\t:";
  if (m_SNPID.size() > 0)
    outfile << '\t' << m_SNPID[m_SNPID.size() - 1];
  if (m_chromosome.size() > 0)
    outfile << '\t' << m_chromosome[m_chromosome.size() - 1];
  if (m_bp.size() > 0)
    outfile << '\t' << m_bp[m_bp.size() - 1];
  if (m_refAllele.size() > 0)
    outfile << '\t' << m_refAllele[m_refAllele.size() - 1];
  if (m_altAllele.size() > 0)
    outfile << '\t' << m_altAllele[m_altAllele.size() - 1];
  outfile << std::endl;
  outfile << "First altFreq\t\t:";
  if (m_altFreq.size() != 0) {
    for (dit = m_altFreq[0].begin(); dit != m_altFreq[0].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "Last altFreq\t\t:";
  if (m_altFreq.size() != 0) {
    for (dit = m_altFreq[m_numSNPs - 1].begin(); dit != m_altFreq[m_numSNPs - 1].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "First maf\t\t:";
  if (m_maf.size() != 0) {
    for (dit = m_maf[0].begin(); dit != m_maf[0].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "Last maf\t\t:";
  if (m_maf.size() != 0) {
    for (dit = m_maf[m_numSNPs - 1].begin(); dit != m_maf[m_numSNPs - 1].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "First avgCall\t\t:";
  if (m_avgCall.size() != 0) {
    for (dit = m_avgCall[0].begin(); dit != m_avgCall[0].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "Last avgCall\t\t:";
  if (m_avgCall.size() != 0) {
    for (dit = m_avgCall[m_numSNPs - 1].begin(); dit != m_avgCall[m_numSNPs - 1].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "First r-squared\t\t:";
  if (m_rSq.size() != 0) {
    for (dit = m_rSq[0].begin(); dit != m_rSq[0].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "Last r-squared\t\t:";
  if (m_rSq.size() != 0) {
    for (dit = m_rSq[m_numSNPs - 1].begin(); dit != m_rSq[m_numSNPs - 1].end(); ++dit)
      outfile << '\t' << *dit;
  }
  outfile << std::endl;
  outfile << "Start Dosages\t\t:\t" << m_startDosages << std::endl;
  outfile << "Current SNP\t\t:\t" << m_currentSNP << std::endl;
  outfile << "First Subject SNP\t:";
  if (m_dosage.size() != 0)
    outfile << '\t' << m_dosage[0];
  if (m_p0.size() != 0)
    outfile << '\t' << m_p0[0];
  if (m_p1.size() != 0)
    outfile << '\t' << m_p1[0];
  if (m_p2.size() != 0)
    outfile << '\t' << m_p2[0];
  outfile << std::endl;
  outfile << "Last Subject SNP\t:";
  if (m_dosage.size() != 0)
    outfile << '\t' << m_dosage[m_numSubjects - 1];
  if (m_p0.size() != 0)
    outfile << '\t' << m_p0[m_numSubjects - 1];
  if (m_p1.size() != 0)
    outfile << '\t' << m_p1[m_numSubjects - 1];
  if (m_p2.size() != 0)
    outfile << '\t' << m_p2[m_numSubjects - 1];
  outfile << std::endl;
  outfile << "---------------------------------------------------------" << std::endl;
}

int CReadBinaryDosageX::GetFirst() {
  if (!m_good || m_numSNPs == 0) {
    m_currentSNP = -1;
    return 1;
  }
  m_infile.seekg(m_startDosages);

  m_currentSNP = 0;
  if (m_subVersion == 1) {
    return ReadDosageSNP();
  } else if (m_mainVersion < 3) {
    return ReadFullGeneticSNP();
  }
  return ReadV4GeneticSNP();
}
///////////////////////////////////////////////////////////////////////////////
//
//                            CReadMultifileBinaryDosage
//
///////////////////////////////////////////////////////////////////////////////

CReadMultifileBinaryDosage::CReadMultifileBinaryDosage(const std::vector<std::string> &filenames) : CReadBinaryDosageX(filenames) {
  if (m_good == true) {
    if (filenames.size() != 3) {
      m_good = false;
    } else {
      m_famFilename = filenames[1];
      m_mapFilename = filenames[2];
      m_famFile.open(m_famFilename.c_str());
      m_mapFile.open(m_mapFilename.c_str());
      if (!m_famFile.good() || !m_mapFile.good())
        m_good = false;
    }
  }
}

CReadMultifileBinaryDosage::~CReadMultifileBinaryDosage() {
  m_famFile.close();
  m_mapFile.close();
}

int CReadMultifileBinaryDosage::ReadSubjects() {
  std::string junk, fam, sub;
  std::istringstream iss;
  bool familyID;
  int numCol;


  if (!m_famFile.good())
    return 1;

  std::getline(m_famFile, junk);
  if (m_famFile.fail())
    return 1;
  numCol = 0;
  iss.str(junk);
  iss >> fam;
  while (!iss.fail()) {
    ++numCol;
    iss >> fam;
  }
  switch(numCol) {
  case 1:
  case 5:
    familyID = false;
    break;
  case 2:
  case 6:
    familyID = true;
    break;
  default:
    return 1;
  }
  do {
    iss.clear();
    iss.str(junk);
    if (familyID) {
      iss >> fam;
      m_FID.push_back(fam);
    }
    iss >> sub;
    m_IID.push_back(sub);
    getline(m_famFile, junk);
  } while (!m_famFile.fail());
  m_numSubjects = m_IID.size();
  return 0;
}

int CReadMultifileBinaryDosage::ReadSNPs() {
  std::string junk, chromosome, snpID, refAllele, altAllele;
  int bp;
  std::istringstream iss;


  if (!m_mapFile.good())
    return 1;

  std::getline(m_mapFile, junk);
  if (m_mapFile.fail())
    return 1;
  iss.str(junk);
  iss >> chromosome >> snpID >> bp >> bp >> refAllele >> altAllele;
  if (iss.fail())
    return 1;
  do {
    iss.clear();
    iss.str(junk);
    iss >> chromosome >> snpID >> bp >> bp >> refAllele >> altAllele;
    m_chromosome.push_back(chromosome);
    m_SNPID.push_back(snpID);
    m_bp.push_back(bp);
    m_refAllele.push_back(refAllele);
    m_altAllele.push_back(altAllele);
    getline(m_mapFile, junk);
  } while (!m_mapFile.fail());
  m_numSNPs = m_SNPID.size();
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage11
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage11::CReadBinaryDosage11(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 1;
  m_subVersion = 1;
  m_scale = 32767.;
}

int CReadBinaryDosage11::ReadHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x01};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage12
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage12::CReadBinaryDosage12(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 1;
  m_subVersion = 2;
  m_scale = 65534.;
}

int CReadBinaryDosage12::ReadHeader() {
  const char version[4] = { 0x00, 0x01, 0x00, 0x02};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage21
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage21::CReadBinaryDosage21(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 2;
  m_subVersion = 1;
}

int CReadBinaryDosage21::ReadHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x01};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage22
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage22::CReadBinaryDosage22(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 2;
  m_subVersion = 2;
}

int CReadBinaryDosage22::ReadHeader() {
  const char version[4] = { 0x00, 0x02, 0x00, 0x02};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage31
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage31::CReadBinaryDosage31(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 3;
  m_subVersion = 1;
  m_startDosages = 12;
}

int CReadBinaryDosage31::ReadHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x01};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

int CReadBinaryDosage31::ReadSubjects() {
  m_infile.read((char *)&m_numSubjects, sizeof(int));
  if (!m_infile.good())
    return 1;
  m_IID.reserve(m_numSubjects);
  m_FID.reserve(m_numSubjects);
  return CReadMultifileBinaryDosage::ReadSubjects();
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage32
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage32::CReadBinaryDosage32(const std::vector<std::string> &filenames) : CReadMultifileBinaryDosage(filenames) {
  m_mainVersion = 3;
  m_subVersion = 2;
  m_startDosages = 12;
}

int CReadBinaryDosage32::ReadHeader() {
  const char version[4] = { 0x00, 0x03, 0x00, 0x02};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  return ReadVersion(version);
}

int CReadBinaryDosage32::ReadSubjects() {
  m_infile.read((char *)&m_numSubjects, sizeof(int));
  if (!m_infile.good())
    return 1;
  m_IID.reserve(m_numSubjects);
  m_FID.reserve(m_numSubjects);
  return CReadMultifileBinaryDosage::ReadSubjects();
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage4x
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage4x::CReadBinaryDosage4x(const std::vector<std::string> &filenames) : CReadBinaryDosageX(filenames) {
  if (filenames.size() != 1)
    m_good = false;
}

int CReadBinaryDosage4x::ReadString(std::vector<std::string> &stringToRead, const int sizeToRead) {
  std::vector<std::string>::iterator stringIt;
  std::string x;
  char *charString = NULL;
  std::istringstream iss;

  charString = new char[sizeToRead];
  m_infile.read(charString, sizeToRead);
  x = charString;
  iss.str(x);
  for (stringIt = stringToRead.begin(); stringIt != stringToRead.end(); ++stringIt) {
    iss >> x;
    *stringIt = x;
  }

  if (charString)
    delete [] charString;

  if (iss.fail())
    return 1;

  return 0;
}

int CReadBinaryDosage4x::ReadHeader() {
  m_infile.seekg((int)Header4pos::numSub);
  m_infile.read((char *)&m_numSubjects, sizeof(int));
  m_infile.read((char *)&m_numSNPs, sizeof(int));
  m_infile.read((char *)&m_numGroups, sizeof(int));
  m_infile.read((char *)&m_subjectOptions, sizeof(int));
  m_infile.read((char *)&m_snpOptions, sizeof(int));
  m_infile.read((char *)&m_startSubjects, sizeof(int));
  m_infile.read((char *)&m_startSNPs, sizeof(int));
  m_infile.read((char *)&m_startDosages, sizeof(int));
  if (m_infile.fail())
    return 1;
  return 0;
}

int CReadBinaryDosage4x::ReadGroups() {
  m_groupSize.resize(m_numGroups);
  m_infile.seekg((int)Header4pos::startGroups);
  m_infile.read((char *)m_groupSize.data(), m_numGroups * sizeof(int));
  if (!m_infile.good())
    return 1;

  return 0;
}

int CReadBinaryDosage4x::ReadSubjects() {
  int iidSize, fidSize;

  m_infile.seekg(m_startSubjects);
  m_infile.read((char *)&iidSize, sizeof(int));
  m_infile.read((char *)&fidSize, sizeof(int));

  m_IID.resize(m_numSubjects);
  if (ReadString(m_IID, iidSize))
    return 1;

  if (fidSize > 0) {
    m_FID.resize(m_numSubjects);
    if (ReadString(m_FID, fidSize))
      return 1;
  } else {
    m_FID.resize(0);
  }

  return 0;
}

int CReadBinaryDosage4x::ReadSNPs() {
  int sizeID, sizeChr, sizeRef, sizeAlt;
  std::vector<std::vector<double> >::iterator dvIt;

  m_infile.seekg(m_startSNPs);
  m_infile.read((char *)&sizeID, sizeof(int));
  m_infile.read((char *)&sizeChr, sizeof(int));
  m_infile.read((char *)&sizeRef, sizeof(int));
  m_infile.read((char *)&sizeAlt, sizeof(int));

  if (m_snpOptions & 0x0002) {
    m_SNPID.resize(m_numSNPs);
    ReadString(m_SNPID, sizeID);
  }
  if (m_snpOptions & 0x0004) {
    if (m_snpOptions & 0x0008)
      m_chromosome.resize(1);
    else
      m_chromosome.resize(m_numSNPs);
    ReadString(m_chromosome, sizeChr);
  }
  m_bp.resize(m_numSNPs);
  m_infile.read((char *)m_bp.data(), m_numSNPs * sizeof(int));
  if (m_snpOptions & 0x0020) {
    m_refAllele.resize(m_numSNPs);
    ReadString(m_refAllele, sizeRef);
  }
  if (m_snpOptions & 0x0040) {
    m_altAllele.resize(m_numSNPs);
    ReadString(m_altAllele, sizeAlt);
  }
  if (m_snpOptions & 0x0080) {
    m_altFreq.resize(m_numSNPs);
    for (dvIt = m_altFreq.begin(); dvIt != m_altFreq.end(); ++dvIt) {
      dvIt->resize(m_numGroups);
      m_infile.read((char *)dvIt->data(), m_numGroups * sizeof(double));
    }
  }
  if (m_snpOptions & 0x0100) {
    m_maf.resize(m_numSNPs);
    for (dvIt = m_maf.begin(); dvIt != m_maf.end(); ++dvIt) {
      dvIt->resize(m_numGroups);
      m_infile.read((char *)dvIt->data(), m_numGroups * sizeof(double));
    }
  }
  if (m_snpOptions & 0x0200) {
    m_avgCall.resize(m_numSNPs);
    for (dvIt = m_avgCall.begin(); dvIt != m_avgCall.end(); ++dvIt) {
      dvIt->resize(m_numGroups);
      m_infile.read((char *)dvIt->data(), m_numGroups * sizeof(double));
    }
  }
  if (m_snpOptions & 0x0400) {
    m_rSq.resize(m_numSNPs);
    for (dvIt = m_rSq.begin(); dvIt != m_rSq.end(); ++dvIt) {
      dvIt->resize(m_numGroups);
      m_infile.read((char *)dvIt->data(), m_numGroups * sizeof(double));
    }
  }
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//                            CReadBinaryDosage41
//
///////////////////////////////////////////////////////////////////////////////

CReadBinaryDosage41::CReadBinaryDosage41(const std::vector<std::string> &filenames) : CReadBinaryDosage4x(filenames) {
  m_mainVersion = 4;
  m_subVersion = 1;
  m_startDosages = 0;
}

int CReadBinaryDosage41::ReadHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x01};

  if (CReadBinaryDosageX::ReadHeader())
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

CReadBinaryDosage42::CReadBinaryDosage42(const std::vector<std::string> &filenames) : CReadBinaryDosage4x(filenames) {
  m_mainVersion = 4;
  m_subVersion = 2;
  m_startDosages = 0;
}

int CReadBinaryDosage42::ReadHeader() {
  const char version[4] = { 0x00, 0x04, 0x00, 0x02};

  if (CReadBinaryDosageX::ReadHeader())
    return 1;
  if (ReadVersion(version))
    return 1;
  return CReadBinaryDosage4x::ReadHeader();
}

///////////////////////////////////////////////////////////////////////////////
//                            Test code
///////////////////////////////////////////////////////////////////////////////

int TestReadBD(CReadBinaryDosageX *bdf) {
  bdf->ReadHeader();
  bdf->ReadSubjects();
  bdf->ReadGroups();
  bdf->ReadSNPs();
  bdf->GetFirst();
  bdf->WriteData(Rcpp::Rcout);
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
  std::vector<std::string> f11Files, f12Files, f21Files, f22Files, f31Files, f32Files, f41Files, f42Files;

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

  CReadBinaryDosage11 f11(f11Files);
  CReadBinaryDosage12 f12(f12Files);
  CReadBinaryDosage21 f21(f21Files);
  CReadBinaryDosage22 f22(f22Files);
  CReadBinaryDosage31 f31(f31Files);
  CReadBinaryDosage32 f32(f32Files);
  CReadBinaryDosage41 f41(f41Files);
  CReadBinaryDosage42 f42(f42Files);

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
