#ifndef BDOSEMINIREADER_H
#define BDOSEMINIREADER_H 1

#ifndef GENETICDATAREADER_H
#include "GeneticDataReader.h"
#endif

int GetBDoseFormat(const std::string &_filename, int &_version, int &_subversion);

class CBDoseMiniReader {
protected:
  bool m_good;
  std::string m_filename;
  std::ifstream m_infile;
  int m_version, m_subversion;
  unsigned int m_numSamples, m_numSNPs;

  std::streampos m_startDosageData;

  CGeneticDataReader *m_geneticDataReader;
  std::vector<double> m_dosage, m_p0, m_p1, m_p2;
  unsigned int m_currentSNP;

  CBDoseMiniReader(const std::string &_filename);
public:
  virtual ~CBDoseMiniReader();

  bool GetFirst();
  bool GetNext();
  bool GetSNP(unsigned int n);
  bool GetSNP(unsigned int n, std::streampos snpLoc);

  bool good() const { return m_good; }
  int Version() const { return m_version;  }
  int SubVersion() const { return m_subversion; }
  unsigned int NumSamples() const { return m_numSamples; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }
  int CurrentSNP() const { return m_currentSNP + 1; }
};

class CBDoseMiniReader1 : public CBDoseMiniReader {
protected:
public:
  CBDoseMiniReader1(const std::string &_filename, const unsigned int _numSamples, const unsigned int _numSNPs);
  virtual ~CBDoseMiniReader1() {};
};

class CBDoseMiniReader4 : public CBDoseMiniReader {
protected:
public:
  CBDoseMiniReader4(const std::string &_filename);
  virtual ~CBDoseMiniReader4() {} ;
};

#endif
