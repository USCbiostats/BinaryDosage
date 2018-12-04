#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "MiniReader.h"

CMiniReader::CMiniReader(const std::string &_filename) {
  m_good = true;
  m_filename = _filename;
  m_numSamples = 0;
  m_numSNPs = 0;
  m_currentPos = 0;
  m_startDosageData = 0;
  m_geneticDataReader = NULL;
  m_currentSNP = 0;
}

CMiniReader::~CMiniReader() {
  if (m_geneticDataReader)
    delete m_geneticDataReader;
  m_infile.close();
}

bool CMiniReader::GetFirst() {
  if (!m_good)
    return false;
  m_currentSNP = 0;
  m_infile.seekg(m_startDosageData);
  if (m_geneticDataReader->ReadData(m_infile, m_dosage, m_p0, m_p1, m_p2))
    m_good = false;
  return m_good;
}

bool CMiniReader::GetNext() {
  if (!m_good)
    return false;
  ++m_currentSNP;
  if (m_currentSNP == NumSNPs()) {
    --m_currentSNP;
    return false;
  }
  if (m_geneticDataReader->ReadData(m_infile, m_dosage, m_p0, m_p1, m_p2))
    m_good = false;
  return m_good;
}

bool CMiniReader::GetSNP(int n) {
  if (!m_good)
    return false;
  if (n == 0 || n > NumSNPs())
    return false;
  --n;
  if (n == m_currentSNP)
    return true;
  if (n == 0)
    return GetFirst();
  if (n < m_currentSNP) {
    m_infile.seekg(m_startDosageData);
    m_geneticDataReader->SkipSNP(m_infile);
    m_currentSNP = 0;
  }
  while (m_currentSNP < n - 1) {
    m_geneticDataReader->SkipSNP(m_infile);
    ++m_currentSNP;
  }
  return GetNext();
}

bool CMiniReader::GetSNP(int n, std::streampos snpLoc) {
  m_infile.seekg(snpLoc);
  m_currentSNP = n - 1;
  if (m_geneticDataReader->ReadData(m_infile, m_dosage, m_p0, m_p1, m_p2))
    m_good = false;
  return m_good;
}

int CMiniReader::CloseFile() {
  if (!m_good)
    return 1;

  if (m_infile.is_open()) {
    m_currentPos = m_infile.tellg();
    m_infile.close();
  }

  return 0;
}
