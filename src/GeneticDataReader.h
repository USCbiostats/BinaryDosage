#ifndef GENETICDATAREADER_H
#define GENETICDATAREADER_H 1

class CGeneticDataReader {
protected:
  const unsigned int m_sampleSize;

  CGeneticDataReader(const unsigned int _sampleSize) : m_sampleSize(_sampleSize) {}
public:
  virtual ~CGeneticDataReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2) = 0;
  virtual int SkipSNP(std::ifstream &_infile) = 0;
};

class CBDoseDataReader : public CGeneticDataReader {
protected:
  const unsigned short m_scale;
  double m_dScale;

  std::vector<unsigned short> m_dataToRead;

  CBDoseDataReader(const unsigned short _scale, const unsigned int _sampleSize);
public:
  virtual ~CBDoseDataReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2) = 0;
  virtual int SkipSNP(std::ifstream &_infile) = 0;
};

class CBDoseDosageReader : public CBDoseDataReader {
protected:
public:
  CBDoseDosageReader(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CBDoseDosageReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

class CBDose1DataReader : public CBDoseDataReader {
protected:
public:
  CBDose1DataReader(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CBDose1DataReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

class CBDose3DataReader : public CBDoseDataReader {
protected:
public:
  CBDose3DataReader(const unsigned short _scale, const unsigned int _sampleSize);
  virtual ~CBDose3DataReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

class CVCFDataReader : public CGeneticDataReader {
protected:
public:
  CVCFDataReader(const unsigned int _sampleSize) : CGeneticDataReader(_sampleSize) {}
  virtual ~CVCFDataReader() {}

  virtual int ReadData(std::ifstream &_infile,
                       std::vector<double> &_dosage,
                       std::vector<double> &_p0,
                       std::vector<double> &_p1,
                       std::vector<double> &_p2);
  virtual int SkipSNP(std::ifstream &_infile);
};

#endif
