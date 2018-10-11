#ifndef GENETICDATAWRITER_H
#define GENETICDATAWRITER_H 1

class CGeneticDataWriter {
protected:
  const unsigned short m_scale;
  double m_dScale;

  const int m_sampleSize;

  std::vector<unsigned short> m_dataToWrite;

  short ConvertToShort(const double x);
  CGeneticDataWriter(const unsigned short _scale, const int _sampleSize);
public:
  virtual ~CGeneticDataWriter() {};

  virtual int WriteData(std::fstream &_outfile,
                        const std::vector<double> &_dosage,
                        const std::vector<double> &_p0,
                        const std::vector<double> &_p1,
                        const std::vector<double> &_p2) = 0;
};

class CDosageDataWriter : public CGeneticDataWriter {
protected:
public:
  CDosageDataWriter(const unsigned short _scale, const int _sampleSize);
  virtual ~CDosageDataWriter() {};

  virtual int WriteData(std::fstream &_outfile,
                        const std::vector<double> &_dosage,
                        const std::vector<double> &_p0,
                        const std::vector<double> &_p1,
                        const std::vector<double> &_p2);
};

class CGeneticDataWriter1 : public CGeneticDataWriter {
protected:
public:
  CGeneticDataWriter1(const unsigned short _scale, const int _sampleSize);
  virtual ~CGeneticDataWriter1() {};

  virtual int WriteData(std::fstream &_outfile,
                        const std::vector<double> &_dosage,
                        const std::vector<double> &_p0,
                        const std::vector<double> &_p1,
                        const std::vector<double> &_p2);
};

class CGeneticDataWriter3 : public CGeneticDataWriter {
protected:
public:
  CGeneticDataWriter3(const unsigned short _scale, const int _sampleSize);
  virtual ~CGeneticDataWriter3() {};

  virtual int WriteData(std::fstream &_outfile,
                        const std::vector<double> &_dosage,
                        const std::vector<double> &_p0,
                        const std::vector<double> &_p1,
                        const std::vector<double> &_p2);
};

#endif
