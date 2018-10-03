#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GeneticDataReader.h"
#include "VCFMiniReader.h"

CVCFMiniReader::CVCFMiniReader(const std::string &_filename, const unsigned int _numSamples,
                               const unsigned int _numSNPs, const std::streampos _startData) : CMiniReader(_filename) {
  m_good = false;
  m_numSamples = _numSamples;
  m_numSNPs = _numSNPs;
  m_startDosageData = _startData;

  m_infile.open(m_filename.c_str());
  if (!m_infile.good())
    return;
  m_good = true;
}
