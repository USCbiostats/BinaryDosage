#ifndef READBINARYDOSAGE_H
#define READBINARYDOSAGE_H 1

class CReadBinaryDosageX {
protected:
  std::string m_filename;
  std::ifstream m_infile;
  char m_version[4];
  int m_mainVersion, m_subVersion;
  int m_numSubjects, m_numSNPs, m_numGroups;
  int m_subjectOptions, m_snpOptions;
  std::vector<int> m_groupSize;
  std::vector<std::string> m_FID, m_IID;
  std::vector<std::string> m_SNPID, m_chromosome, m_refAllele, m_altAllele;
  std::vector<int> m_bp;
  bool m_usesFamilyID;

  int ReadVersion();
  CReadBinaryDosageX(const std::string &filename);
public:
  virtual ~CReadBinaryDosageX();

  virtual int ReadHeader();
  virtual int ReadGroups() { return 0; }
  virtual int ReadSubjects() = 0;
  virtual int ReadSNP() = 0;

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
};

class CReadMultifileBinaryDosage : public CReadBinaryDosageX {
protected:
  std::ifstream m_famFile;
  std::ifstream m_mapFile;

  CReadMultifileBinaryDosage(const std::string &filename);
public:
  virtual ~CReadMultifileBinaryDosage();

  virtual int ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int ReadSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele);
};

class CReadBinaryDosage11 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage11(const std::string &filename);
  virtual ~CReadBinaryDosage11() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage12 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage12(const std::string &filename);
  virtual ~CReadBinaryDosage12() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage21 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage21(const std::string &filename);
  virtual ~CReadBinaryDosage21() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage22 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage22(const std::string &filename);
  virtual ~CReadBinaryDosage22() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage31 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage31(const std::string &filename);
  virtual ~CReadBinaryDosage31() {}

  virtual int ReadHeader();
  virtual int ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
};

class CReadBinaryDosage32 : public CReadMultifileBinaryDosage {
public:
  CReadBinaryDosage32(const std::string &filename);
  virtual ~CReadBinaryDosage32() {}

  virtual int ReadHeader();
  virtual int ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
};

class CReadBinaryDosage4x : public CReadBinaryDosageX {
protected:
  std::vector<std::string> m_chromosome, m_snpID, m_refAllele, m_altAllele;
  std::vector<int> m_bp;
  std::vector<double> m_calculatedMAF;
  CReadBinaryDosage4x(const std::string &filename);
public:
  virtual ~CReadBinaryDosage4x() {}

  virtual int ReadHeader();
  virtual int ReadGroups(const std::vector<int> &groupSize);
  virtual int ReadSubjects(const std::vector<std::string> &FID, const std::vector<std::string> &IID);
  virtual int ReadSNP(const std::string &chromosome, const std::string &snpID, const int bp,
                       const std::string &refAllele, const std::string &altAllele);
};

class CReadBinaryDosage41 : public CReadBinaryDosage4x {
public:
  CReadBinaryDosage41(const std::string &filename);
  virtual ~CReadBinaryDosage41() {}

  virtual int ReadHeader();
};

class CReadBinaryDosage42 : public CReadBinaryDosage4x {
public:
  CReadBinaryDosage42(const std::string &filename);
  virtual ~CReadBinaryDosage42() {}

  virtual int ReadHeader();
};

#endif
