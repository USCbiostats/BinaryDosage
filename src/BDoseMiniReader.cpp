#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GeneticDataReader.h"
#include "BDoseMiniReader.h"

//#define BINARYDOSAGEMINIREADER_DEBUG 1

const char BDHeader[4] = { 'b', 'o', 's', 'e' };
const int BDHeaderSize = 8;

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

CBDoseMiniReader::CBDoseMiniReader(const std::string &_filename) : CMiniReader(_filename) {
  m_good = false;
  m_version = 0;
  m_subversion = 0;
  if (GetBDoseFormat(m_filename, m_version, m_subversion))
    return;
  m_infile.open(_filename.c_str(), std::ios_base::in | std::ios_base::binary);
  if (!m_infile.good())
    return;
  m_good = true;
}


CBDoseMiniReader1::CBDoseMiniReader1(const std::string &_filename, const unsigned int _numSamples, const unsigned int _numSNPs) : CBDoseMiniReader(_filename) {
  if (!m_good)
    return;

  if (m_version == 3) {
    m_infile.read((char *)&m_numSamples, sizeof(int));
    if (m_numSamples == _numSamples)
      return;
    m_startDosageData = 12;
  } else {
    m_numSamples = _numSamples;
    m_startDosageData = 8;
  }
  m_numSNPs = _numSNPs;
  m_good = true;
}

CBDoseMiniReader4::CBDoseMiniReader4(const std::string &_filename) : CBDoseMiniReader(_filename) {
  int startDosageData;

  if (!m_good)
    return;
  m_good = false;

  if (m_version != 4)
    return;

  m_infile.seekg((int)Header4pos::numSub);
  m_infile.read((char *)&m_numSamples, sizeof(int));
  m_infile.read((char *)&m_numSNPs, sizeof(int));
  m_infile.seekg((int)Header4pos::startDosage);
  m_infile.read((char *)&startDosageData, sizeof(int));
#ifdef BINARYDOSAGEMINIREADER_DEBUG
  std::cout << m_numSamples << '\t' << m_numSNPs << '\t' << startDosageData << std::endl;
# endif
  m_startDosageData = startDosageData;

  if (!m_infile.good())
    return;

  m_dosage.resize(NumSamples());
  m_p0.resize(NumSamples());
  m_p1.resize(NumSamples());
  m_p2.resize(NumSamples());

  if (m_subversion == 1)
    m_geneticDataReader = new CBDoseDosageReader(10000, NumSamples());
  else
    m_geneticDataReader = new CBDose3DataReader(10000, NumSamples());

  m_good = true;
  m_good = GetFirst();
}
