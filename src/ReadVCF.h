#ifndef READVCF_H
#define READVCF_H 1

class CReadVCFBase {
protected:
  std::string m_filename;
  std::ifstream m_infile;
  
  int m_subjectLine;
  int m_startLine;
  
  int m_numSubjects;
  std::vector<std::string> m_subjectID;
  
  int m_numSNPs;
  std::vector<std::string> m_snpID;
  std::vector<std::string> m_chromosome;
  std::vector<int> m_bp;
  std::vector<std::string> m_refAllele;
  std::vector<std::string> m_altAllele;
  std::vector<double> m_altAlleleFreq;
  std::vector<double> m_minorAlleleFreq;
  std::vector<double> m_averageCall;
  std::vector<double> m_rSquared;
  
  std::vector<double> m_dosage;
  std::vector<double> m_p0;
  std::vector<double> m_p1;
  std::vector<double> m_p2;
  
  int GetFirst();
  int GetNext();

public:
  CReadVCFBase(const std::string &filename);
  virtual ~CReadVCFBase();

  virtual int ReadSubjects();
  virtual int GetFirstSNP() = 0;
  virtual int GetNextSNP() = 0;
  
  double AlternateAlleleFrequency();

  int SubjectLine() const { return m_subjectLine; }
  int StartLine() const { return m_startLine; }
  
  int NumSubjects() const { return m_numSubjects; }
  const std::vector<std::string> &SubjectID() const { return m_subjectID; }
  
  int NumSNPs() const { return m_numSNPs; }
  const std::vector<std::string> &SNPID() const { return m_snpID; }
  const std::vector<std::string> &Chromosome() const { return m_chromosome; }
  const std::vector<int> &Location() const { return m_bp; }
  const std::vector<std::string> &ReferenceAllele() const { return m_refAllele; }
  const std::vector<std::string> &AlternateAllele() const { return m_altAllele; }
  const std::vector<double> &AltAlleleFreq() const { return m_altAlleleFreq; }
  const std::vector<double> &MAF() const { return m_minorAlleleFreq; }
  const std::vector<double> &AverageCall() const { return m_averageCall; }
  const std::vector<double> &RSquared() const { return m_rSquared; }
  
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }
};

class CReadVCF_HRC : public CReadVCFBase {
protected:
  double ReadDosage(const char *x);
  double ReadProbability(const char *x);
public:
  CReadVCF_HRC(const std::string &filename) : CReadVCFBase(filename) {}
  virtual ~CReadVCF_HRC() {}
  
  virtual int GetFirstSNP();
  virtual int GetNextSNP();
};
#endif