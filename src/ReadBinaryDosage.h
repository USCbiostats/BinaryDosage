#ifndef READBINARYDOSAGE_H
#define READBINARYDOSAGE_H 1

class CReadBinaryDosageX {
protected:
  std::string m_filename;
  std::ifstream m_infile;
  bool m_good;
  char m_version[4];
  int m_mainVersion, m_subVersion;
  int m_numSubjects, m_numSNPs, m_numGroups;
  int m_startSubjects, m_startSNPs, m_startDosages;
  int m_subjectOptions, m_snpOptions;
  int m_currentSNP;
  std::vector<int> m_groupSize;
  std::vector<std::string> m_FID, m_IID;
  std::vector<std::string> m_SNPID, m_chromosome, m_refAllele, m_altAllele;
  std::vector<int> m_bp;
  std::vector<std::vector<double> > m_altFreq, m_maf, m_avgCall, m_rSq;
  bool m_usesFamilyID;
  double m_scale;
  std::vector<unsigned short> m_sReadValues;
  std::vector<double> m_dosage, m_p0, m_p1, m_p2;

  int ReadVersion(const char *version);
//  virtual int ProcessSNP() = 0;
  int ReadDosageSNP();
  int ReadFullGeneticSNP();
  int ReadV4GeneticSNP();
  CReadBinaryDosageX(const std::vector<std::string> &filenames);
public:
  virtual ~CReadBinaryDosageX();

  virtual int ReadHeader();
  virtual int ReadGroups();
  virtual int ReadSubjects() = 0;
  virtual int ReadSNPs() = 0;
  int GetFirst();
  virtual int GetNext();
  virtual int GetSNP(const int n);

  bool good() const { return m_good; }
  int Version() const { return m_mainVersion; }
  int SubVersion() const { return m_subVersion; }
  int NumSubjects() const { return m_numSubjects; }
  int NumSNPs() const { return m_numSNPs; }
  int NumGroups() const { return m_numGroups; }
  const std::vector<int> &GroupSize() const { return m_groupSize; }
  const std::vector<std::string> &FamilyID() const { return m_FID; }
  const std::vector<std::string> &SubjectID() const { return m_IID; }
  const std::vector<std::string> &SNPID() const { return m_SNPID; }
  const std::vector<std::string> &Chromosome() const { return m_chromosome; }
  const std::vector<int> &Location() const { return m_bp; }
  const std::vector<std::string> &ReferenceAllele() const { return m_refAllele; }
  const std::vector<std::string> &AlternateAllele() const { return m_altAllele; }
  const std::vector<std::vector<double> > &AltAlleleFreq() const { return m_altFreq; }
  const std::vector<std::vector<double> > &MAF() const { return m_maf; }
  const std::vector<std::vector<double> > &AvgCall() const { return m_avgCall; }
  const std::vector<std::vector<double> > &RSquared() const { return m_rSq; }
  const std::vector<unsigned short> &ReadValues() const { return m_sReadValues; }
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }

  void WriteData(std::ostream &outfile);
};

class CReadMultifileBinaryDosage : public CReadBinaryDosageX {
protected:
  std::string m_famFilename;
  std::string m_mapFilename;
  std::ifstream m_famFile;
  std::ifstream m_mapFile;

  CReadMultifileBinaryDosage(const std::vector<std::string> &filenames);
public:
  virtual ~CReadMultifileBinaryDosage();

  virtual int ReadSubjects();
  virtual int ReadSNPs();
};

class CReadBinaryDosage11 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage11(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage11() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage12 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage12(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage12() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage21 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage21(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage21() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage22 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage22(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage22() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage31 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage31(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage31() {}

  virtual int ReadHeader();
  virtual int ReadSubjects();
};

class CReadBinaryDosage32 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage32(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage32() {}

  virtual int ReadHeader();
  virtual int ReadSubjects();
};

class CReadBinaryDosage4x : public CReadBinaryDosageX {
protected:
  int ReadString(std::vector<std::string> &stringToRead, const int sizetoRead);
  CReadBinaryDosage4x(const std::vector<std::string> &filenames);
public:
  virtual ~CReadBinaryDosage4x() {}

  virtual int ReadHeader();
  virtual int ReadGroups();
  virtual int ReadSubjects();
  virtual int ReadSNPs();
};

class CReadBinaryDosage41 : public CReadBinaryDosage4x {
public:
  CReadBinaryDosage41(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage41() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage42 : public CReadBinaryDosage4x {
public:
  CReadBinaryDosage42(const std::vector<std::string> &filenames);
  virtual ~CReadBinaryDosage42() {}

  virtual int ReadHeader();
};

#endif
