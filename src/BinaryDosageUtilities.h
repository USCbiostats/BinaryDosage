#ifndef BINARYDOSAGEUTILITIES_H
#define BINARYDOSAGEuTILITIES_H 1

class CReadBinaryDosage {
protected:
  static const char m_magicWord42[8];
  char m_magicWord[8];
  unsigned int m_numSubjects;
  unsigned int m_numSNPs;
  unsigned int m_numGroups;
  unsigned int m_subjectOptions;
  unsigned int m_SNPOptions;
  unsigned int m_subjectStart;
  unsigned int m_SNPStart;
  unsigned int m_dosageStart;
  std::vector <unsigned int> m_groupSize;

  std::vector<std::string> m_SNPName;
  std::vector<std::string> m_chromosome;
  std::vector<unsigned int> m_location;
  std::vector<std::string> m_refAllele;
  std::vector<std::string> m_altAllele;
  std::vector<double> m_altFreq;
  std::vector<double> m_maf;
  std::vector<double> m_avgCall;
  std::vector<double> m_rSquared;
  unsigned int m_SNPNameSize;
  unsigned int m_chromosomeSize;
  unsigned int m_refAlleleSize;
  unsigned int m_altAlleleSize;

  std::vector <unsigned int> m_SNPDataSize;

  std::ifstream m_infile;

  std::vector<std::string> m_subjectID;
  std::vector<std::string> m_familyID;
  unsigned int m_subjectIDSize;
  unsigned int m_familyIDSize;

  int ReadHeader();
  int ReadSubjects();
  int ReadSNPInfo();
  int ReadSNPDataSize();

public:
  CReadBinaryDosage();
  ~CReadBinaryDosage();

  int ReadFileInfo(const std::string &binaryDosageFilename);
//  int ReadDosage(unsigned int n, Rcpp::NumericVector &d, unsigned int numSteps, unsigned int dosageStart);

  unsigned int NumSubjects() const { return m_numSubjects; }
  unsigned int NumSNPs() const { return m_numSNPs; }
  unsigned int NumGroups() const { return m_numGroups; }
  unsigned int SubjectOptions() const { return m_subjectOptions; }
  unsigned int SNPOptions() const { return m_SNPOptions; }
  const std::vector<unsigned int> &GroupSize() const { return m_groupSize; }
  unsigned int SubjectStart() const { return m_subjectStart; }
  unsigned int SNPStart() const { return m_SNPStart; }
  unsigned int DosageStart() const { return m_dosageStart; }
  const std::vector<std::string> &SubjectID() { return m_subjectID; }
  const std::vector<std::string> &FamilyID() { return m_familyID; }
  const std::vector<std::string> &SNPNames() { return m_SNPName; }
  const std::vector<std::string> &Chromosome() { return m_chromosome; }
  const std::vector<unsigned int> &Location() { return m_location; }
  const std::vector<std::string> &RefAllele() { return m_refAllele; }
  const std::vector<std::string> &AltAllele() { return m_altAllele; }
  const std::vector<double> &AltFreq() { return m_altFreq; }
  const std::vector<double> &MAF() { return m_maf; }
  const std::vector<double> &AvgCall() { return m_avgCall; }
  const std::vector<double> &RSquared() { return m_rSquared; }
  const std::vector<unsigned int> &SNPDataSize() { return m_SNPDataSize; }
};
#endif
