#ifndef BINARYDOSAGEMINIREADER_H
#define BINARYDOSAGEMINIREADER_H 1

#ifndef GENETICDATAREADER_H
#include "GeneticDataReader.h"
#endif

int GetBinaryDosageFormat(const std::string &_filename, int &_version, int &_subversion);

int ReadFamFile(const std::string &_filename, std::vector<std::string> &_FID, std::vector<std::string> &_SID);
int ReadMapFile(const std::string &_filename, std::vector<std::string> &_chromosome, std::vector<std::string> &_snpID,
                std::vector<int> &_location, std::vector<std::string> &_refAllele, std::vector<std::string> &_altAllele);

class CBinaryDosageMiniReader {
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

  CBinaryDosageMiniReader(const std::string &_filename);
public:
  virtual ~CBinaryDosageMiniReader();

  bool GetFirst();
  bool GetNext();
  bool GetSNP(unsigned int n);

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

class CBinaryDosageMiniReader1 : public CBinaryDosageMiniReader {
protected:
public:
  CBinaryDosageMiniReader1(const std::string &_filename, const std::string &_famFilename, const std::string &_mapFilename) : CBinaryDosageMiniReader(_filename) {}
  virtual ~CBinaryDosageMiniReader1() {};
};

class CBinaryDosageMiniReader4 : public CBinaryDosageMiniReader {
protected:
public:
  CBinaryDosageMiniReader4(const std::string &_filename);
  virtual ~CBinaryDosageMiniReader4() {} ;
};

#endif
