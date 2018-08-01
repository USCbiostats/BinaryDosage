#ifndef WRITEBINARYDOSAGE_H
#define WRITEBINARYDOSAGE_H 1

class CWriteBinaryDosage {
protected:
  std::string m_filename;
  std::ofstream m_outfile;

  int WriteVersion(const char *version);
  CWriteBinaryDosage(const std::string &filename);
public:
  virtual ~CWriteBinaryDosage();

  virtual int WriteHeader();
  virtual int WriteGroups(const std::vector<int> &groupSize) { return 0; }
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID) = 0;
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq) = 0;
  virtual int FinalizeSNPs() { return 0; }
};

class CWriteMultifileBinaryDosage : public CWriteBinaryDosage {
protected:
  std::ofstream m_famFile;
  std::ofstream m_mapFile;

  CWriteMultifileBinaryDosage(const std::string &filename);
public:
  virtual ~CWriteMultifileBinaryDosage();

  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq);
};

class CWriteBinaryDosage11 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage11(const std::string &filename);
  virtual ~CWriteBinaryDosage11() {}

  virtual int WriteHeader();
};

class CWriteBinaryDosage12 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage12(const std::string &filename);
  virtual ~CWriteBinaryDosage12() {}

  virtual int WriteHeader();
};

class CWriteBinaryDosage21 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage21(const std::string &filename);
  virtual ~CWriteBinaryDosage21() {}

  virtual int WriteHeader();
};

class CWriteBinaryDosage22 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage22(const std::string &filename);
  virtual ~CWriteBinaryDosage22() {}

  virtual int WriteHeader();
};

class CWriteBinaryDosage31 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage31(const std::string &filename);
  virtual ~CWriteBinaryDosage31() {}

  virtual int WriteHeader();
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
};

class CWriteBinaryDosage32 : public CWriteMultifileBinaryDosage {
public:
  CWriteBinaryDosage32(const std::string &filename);
  virtual ~CWriteBinaryDosage32() {}

  virtual int WriteHeader();
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
};

class CWriteBinaryDosage4x : public CWriteBinaryDosage {
protected:
  std::vector<std::string> m_chromosome, m_snpID, m_refAllele, m_altAllele;
  std::vector<int> m_bp;
  std::vector<double> m_calculatedMAF;
  std::vector<std::vector<double> > m_altFreq, m_maf, m_avgCall, m_rSq;
  int m_numGroups;
  int m_startSubjects, m_startSNPs, m_startDosages;

  int WriteString(const std::vector<std::string> &stringToWrite);
  int AddToStringVector(std::vector<std::string> &addToVector, const std::string &stringToAdd);
  int AddToDoubleVector(std::vector<std::vector<double> > &addToVector, const std::vector<double> &vectorToAdd);
  int GetSNPOptions();
  CWriteBinaryDosage4x(const std::string &filename);
public:
  virtual ~CWriteBinaryDosage4x() {}

  virtual int WriteHeader();
  virtual int WriteGroups(const std::vector<int> &groupSize);
  virtual int WriteSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int AddSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele,
                       const std::vector<double> &altFreq, const std::vector<double> &maf,
                       const std::vector<double> &avgCall, const std::vector<double> &rSq);
  virtual int WriteSNPs();
};

class CWriteBinaryDosage41 : public CWriteBinaryDosage4x {
public:
  CWriteBinaryDosage41(const std::string &filename);
  virtual ~CWriteBinaryDosage41() {}

  virtual int WriteHeader();
};

class CWriteBinaryDosage42 : public CWriteBinaryDosage4x {
public:
  CWriteBinaryDosage42(const std::string &filename);
  virtual ~CWriteBinaryDosage42() {}

  virtual int WriteHeader();
};

#endif
