#ifndef WRITEBINARYDOSAGE_H
#define WRITEBINARYDOSAGE_H 1

class CWriteBinaryDosage {
protected:
  bool m_ready;

  char m_version[4];

  std::string m_filename;
  std::ofstream m_outfile;
  std::vector<short> m_dataToWrite;

  short ConvertToShort(const double x, const double scale);

  int AddDosagesOnly(const std::vector<std::vector<double> > &dosageValues, double scale);
  int AddGeneticValues32(const std::vector<std::vector<double> > &dosageValues);

  CWriteBinaryDosage(const std::vector<std::string> &filenames);
public:
  virtual ~CWriteBinaryDosage();

  virtual int WriteHeader();
  virtual int WriteGroups(const std::vector<int> &groupSize) { return 0; }
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) = 0;
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq) = 0;
  virtual int WriteSNPs() { return 0; }
  virtual int WriteAllSNPs(const std::vector<std::string> &chromosome,
                           const std::vector<std::string> &snpID,
                           const std::vector<int> bp,
                           const std::vector<std::string> &refAllele,
                           const std::vector<std::string> &altAllele,
                           const std::vector<std::vector<double> > &altFreq,
                           const std::vector<std::vector<double> > &maf,
                           const std::vector<std::vector<double> > &avgCall,
                           const std::vector<std::vector<double> > &rSq) = 0;
  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues) = 0;
};

class CWriteMultifileBinaryDosage : public CWriteBinaryDosage {
protected:
  std::string m_famFilename;
  std::string m_mapFilename;
  std::ofstream m_famFile;
  std::ofstream m_mapFile;

  CWriteMultifileBinaryDosage(const std::vector<std::string> &filenames);

  int AddGeneticValues1or2(const std::vector<std::vector<double> > &geneticValues, double scale);
public:
  virtual ~CWriteMultifileBinaryDosage();

  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq);
  virtual int WriteAllSNPs(const std::vector<std::string> &chromosome,
                           const std::vector<std::string> &snpID,
                           const std::vector<int> bp,
                           const std::vector<std::string> &refAllele,
                           const std::vector<std::string> &altAllele,
                           const std::vector<std::vector<double> > &altFreq,
                           const std::vector<std::vector<double> > &maf,
                           const std::vector<std::vector<double> > &avgCall,
                           const std::vector<std::vector<double> > &rSq);
};

class CWriteBinaryDosage11 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage11(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage11() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage12 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage12(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage12() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage21 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage21(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage21() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage22 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage22(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage22() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage31 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage31(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage31() {}

  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage32 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage32(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage32() {}

  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage4x : public CWriteBinaryDosage {
protected:
  bool m_usesTempFiles;
  std::string m_chromosomeFilename;
  std::string m_SNPFilename;
  std::string m_refFilename;
  std::string m_altFilename;
  std::string m_bpFilename;
  std::string m_altFreqFilename;
  std::string m_mafFilename;
  std::string m_avgCallFilename;
  std::string m_rSqFilename;
  std::fstream m_chromosomeFile;
  std::fstream m_SNPFile;
  std::fstream m_refFile;
  std::fstream m_altFile;
  std::fstream m_bpFile;
  std::fstream m_altFreqFile;
  std::fstream m_mafFile;
  std::fstream m_avgCallFile;
  std::fstream m_rSqFile;
  std::vector<std::string> m_chromosome, m_snpID, m_refAllele, m_altAllele;
  std::vector<int> m_bp;
  std::vector<double> m_calculatedMAF;
  std::vector<std::vector<double> > m_altFreq, m_maf, m_avgCall, m_rSq;
  int m_numGroups;
  int m_startSubjects, m_startSNPs, m_startDosages;

  int WriteString(const std::vector<std::string> &stringToWrite);
  int WriteStringToTempFile(const std::string &stringToWrite, std::fstream &outfile);
  int WriteIntToTempFile(const int valueToWrite, std::fstream &outfile);
  int WriteDoubleVectorToTempFile(const std::vector<double> &vectorToWrite, std::fstream &outfile);
  int AddToStringVector(std::vector<std::string> &addToVector, const std::string &stringToAdd);
  int AddToDoubleVector(std::vector<std::vector<double> > &addToVector, const std::vector<double> &vectorToAdd);
  int GetSNPOptions();
  int WriteStringVectorToFile(const std::vector<std::string> &stringToWrite, int sizeLocation);

  CWriteBinaryDosage4x(const std::vector<std::string> &filenames);
public:
  virtual ~CWriteBinaryDosage4x();

  virtual int WriteHeader();
  virtual int WriteGroups(const std::vector<int> &groupSize);
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq);
  virtual int WriteSNPs();
  virtual int WriteAllSNPs(const std::vector<std::string> &chromosome,
                           const std::vector<std::string> &snpID,
                           const std::vector<int> bp,
                           const std::vector<std::string> &refAllele,
                           const std::vector<std::string> &altAllele,
                           const std::vector<std::vector<double> > &altFreq,
                           const std::vector<std::vector<double> > &maf,
                           const std::vector<std::vector<double> > &avgCall,
                           const std::vector<std::vector<double> > &rSq);
};

class CWriteBinaryDosage41 : public CWriteBinaryDosage4x {
public:
  CWriteBinaryDosage41(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage41() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

class CWriteBinaryDosage42 : public CWriteBinaryDosage4x {
public:
  CWriteBinaryDosage42(const std::vector<std::string> &filenames);
  virtual ~CWriteBinaryDosage42() {}

  virtual int AddGeneticValues(const std::vector<std::vector<double> > &geneticValues);
};

#endif
