#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GeneticDataReader.h"
#include "VCFMiniReader.h"

CVCFMiniReader::CVCFMiniReader(const std::string &_filename, const int _numSamples,
                               const int _numSNPs, const std::streampos _startData) : CMiniReader(_filename) {
  m_good = false;
  m_numSamples = _numSamples;
  m_numSNPs = _numSNPs;
  m_startDosageData = _startData;

  m_infile.open(m_filename.c_str());
  if (!m_infile.good())
    return;

  m_dosage.resize(NumSamples());
  m_p0.resize(NumSamples());
  m_p1.resize(NumSamples());
  m_p2.resize(NumSamples());
  m_geneticDataReader = new CVCFDataReader(m_numSamples);
  m_good = true;
}

int CVCFMiniReader::OpenFile() {
  if (!m_good)
    return 1;

  if (m_infile.is_open())
    return 0;

  m_infile.open(m_filename.c_str());
  m_infile.seekg(m_currentPos);
  if (!m_infile.good()) {
    m_good = false;
    return 1;
  }

  return 0;
}
