#ifndef BDOSEWRITER_H
#define BDOSEWRITER_H 1

#ifndef GENETICDATAWRITER_H
#include "GeneticDataWriter.h"
#endif

class CBDoseWriter {
protected:
  bool m_good;

  std::string m_filename;
  std::fstream m_outfile;

  const int m_format, m_version;

  bool m_subjectDataWritten;
  int m_numSubjects;

  bool m_SNPDataWritten;
  int m_numSNPs;
  int m_snpsWritten;

  std::streampos m_startDosageData;

  CGeneticDataWriter *m_geneticDataWriter;

  virtual int OpenGeneticDataWriter();

  CBDoseWriter(const std::string &_filename, const int _format, const int _version);
public:
  virtual ~CBDoseWriter();

  virtual int WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID);
  virtual int WriteSNPData(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                           const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                           const std::vector<std::string> &_altAllele, int _snpOptions);
  virtual int WriteAdditionalSNPData(const std::vector<std::vector<double> > &_aaf,
                                     const std::vector<std::vector<double> > &_maf,
                                     const std::vector<std::vector<double> > &_avgCall,
                                     const std::vector<std::vector<double> > &_rSq) { return 0; }
  int WriteGeneticData(const std::vector<double> &_dosage, const std::vector<double> &_p0,
                       const std::vector<double> &_p1, const std::vector<double> &_p2);
  int Finalize();

  bool good() const { return m_good; }
  int Format() const { return m_format; }
  int Version() const { return m_version; }
};

class CBDoseWriter1 : public CBDoseWriter {
protected:
  std::string m_famFilename;
  std::string m_mapFilename;

  virtual int OpenGeneticDataWriter();
public:
  CBDoseWriter1(const std::string &_filename, const std::string &_famFilename, const std::string &_mapFilename,
                const int _format, const int _version);
  virtual ~CBDoseWriter1() {}

  virtual int WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID);
  virtual int WriteSNPData(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                           const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                           const std::vector<std::string> &_altAllele, int _snpOptions);
};

class CBDoseWriter4 : public CBDoseWriter {
protected:
  int ProcessString(const std::vector<std::string> &_stringsToWrite, unsigned int &_strLength);
  int ProcessExtraSNPData(const std::vector<std::vector<double> > &_dataToWrite, const unsigned int _groupSize);

  unsigned int GetSNPOptions(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                             const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
                             const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
                             const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq);
  virtual int OpenGeneticDataWriter() = 0;
public:
  CBDoseWriter4(const std::string &_filename, const int _format, const int _version);
  virtual ~CBDoseWriter4() {}

  virtual int WriteSubjectData(const std::vector<std::string> &_FID, const std::vector<std::string> &_SID);
  virtual int WriteSNPData(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
                           const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                           const std::vector<std::string> &_altAllele, int _snpOptions);
  virtual int WriteAdditionalSNPData(const std::vector<std::vector<double> > &_aaf,
                                     const std::vector<std::vector<double> > &_maf,
                                     const std::vector<std::vector<double> > &_avgCall,
                                     const std::vector<std::vector<double> > &_rSq) { return 0; }
};

#endif
