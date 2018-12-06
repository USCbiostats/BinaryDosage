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
  int m_numSamples, m_numSNPs;

  bool m_chunked;
  int m_currentChunk;
  std::string m_readBuffer;
//  std::vector<char> m_readBuffer;
  std::vector<int> m_startSNP;
  std::vector<int> m_snpChunk;
  std::vector<std::streampos> m_filePos;
  std::vector<std::streampos> m_stringPos;
//  std::string m_readString;
  std::istringstream m_iss;

  std::streampos m_currentPos;
  std::streampos m_startDosageData;

  CGeneticDataReader *m_geneticDataReader;
  std::vector<double> m_dosage, m_p0, m_p1, m_p2;
  int m_currentSNP;

  virtual int ReadChunk(int n);
  bool GetChunkSNP();
  CMiniReader(const std::string &_filename);
public:
  virtual ~CMiniReader();

  virtual bool GetFirst();
  virtual bool GetNext();
  virtual bool GetSNP(int n);
  virtual bool GetSNP(int n, std::streampos snpLoc);

  virtual int CloseFile();
  virtual int OpenFile() = 0;

  void ChunkIt(const std::vector<std::streampos> &indices);

  bool good() const { return m_good; }
  int NumSamples() const { return m_numSamples; }
  int NumSNPs() const { return m_numSNPs; }
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }
  int CurrentSNP() const { return m_currentSNP + 1; }
};
#endif
