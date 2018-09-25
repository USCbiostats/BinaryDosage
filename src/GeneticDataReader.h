#ifndef GENETICDATAREADER_H
#define GENETICDATAREADER_H 1

class CGeneticDataReader {
protected:
  const unsigned short m_scale;
  double m_dScale;

  const unsigned int m_sampleSize;

  std::vector<unsigned short> m_dataToRead;

  CGeneticDataReader(const unsigned short _scale, const unsigned int _sampleSize);
public:
  virtual ~CGeneticDataReader() {};

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2) = 0;
  virtual int SkipSNP(std::ifstream &_infile) = 0;
};

class CDosageDataReader : public CGeneticDataReader {
protected:
public:
  CDosageDataReader(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CDosageDataReader() {};

  virtual int ReadData(std::fstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

class CGeneticDataReader1 : public CGeneticDataReader {
protected:
public:
  CGeneticDataReader1(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CGeneticDataReader1() {};

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

class CGeneticDataReader3 : public CGeneticDataReader {
protected:
public:
  CGeneticDataReader3(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CGeneticDataReader3() {};

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

#endif
