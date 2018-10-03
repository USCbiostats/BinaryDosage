#ifndef MINIREADER_H
#define MINIREADER_H 1

#ifndef GENETICDATAREADER_H
#include "GeneticDataReader.h"
#endif

class CMiniReader {
protected:
  bool m_good;
  std::string m_filename;
  std::ifstream m_infile;
  unsigned int m_numSamples, m_numSNPs;

  std::streampos m_startDosageData;

  CGeneticDataReader *m_geneticDataReader;
  std::vector<double> m_dosage, m_p0, m_p1, m_p2;
  unsigned int m_currentSNP;

  CMiniReader(const std::string &_filename);
public:
  virtual ~CMiniReader();

  virtual bool GetFirst();
  virtual bool GetNext();
  virtual bool GetSNP(unsigned int n);
  virtual bool GetSNP(unsigned int n, std::streampos snpLoc);

  bool good() const { return m_good; }
  unsigned int NumSamples() const { return m_numSamples; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }
  int CurrentSNP() const { return m_currentSNP + 1; }
};
#endif
