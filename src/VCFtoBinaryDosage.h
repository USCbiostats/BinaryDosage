#ifndef VCFtoBinaryDosage_H
#define VCFtoBinaryDosage_H

#ifndef READVCF_H
#include "ReadVCF.h"
#endif

class CVCFtoBinaryDosage42 {
protected:
  CReadVCFBase *m_vcfFile;

  unsigned short *m_bd;
  int m_numSub;
  int m_numSNPs;
  std::string m_bdFilename;
  std::ofstream m_chromosomeFile;
  std::ofstream m_bpFile;
  std::ofstream m_snpIDFile;
  std::ofstream m_refFile;
  std::ofstream m_altFile;
  std::ofstream m_doseFile;
  std::ofstream m_altFreqFile;
  std::ofstream m_binaryDosageFile;

  void OpenTempFiles();
  void CloseTempFiles();
  virtual void WriteDosage(const std::vector<double> &d, const std::vector<double> &p0,
                   const std::vector<double> &p1, const std::vector<double> &p2);

  void WriteHeader();
  void WriteGroups();
  void WriteSubjects(const std::vector<std::string> &subID);
  void WriteSNPInfo(const int snpOptions, const std::string &singleChromosome);
  void WriteGeneticValues();

  virtual int OpenVCF(const std::string &vcfFilename);
public:
  CVCFtoBinaryDosage42();
  ~CVCFtoBinaryDosage42();

  int Convert(const std::string &vcfFilename, const std::string &bdFilename);
};

class CVCF53toBinaryDosage42 : public CVCFtoBinaryDosage42{
protected:
  virtual void WriteDosage(const std::vector<double> &d, const std::vector<double> &p0,
                           const std::vector<double> &p1, const std::vector<double> &p2);
  virtual int OpenVCF(const std::string &vcfFilename);
public:
  CVCF53toBinaryDosage42() : CVCFtoBinaryDosage42() {}
  ~CVCF53toBinaryDosage42() {}
};

#endif
