#ifndef BDOSEREADER_H
#define BDOSEREADER_H 1

#ifndef GENETICDATAREADER_H
#include "GeneticDataReader.h"
#endif

int GetBDoseFormat(const std::string &_filename, int &_version, int &_subversion);

class CBDoseReader {
protected:
  bool m_good;
  std::string m_filename;
  std::ifstream m_infile;
  int m_format, m_version;

  std::streampos m_startDosageData;
  std::vector<int> m_snpIndex;

  std::vector<int> m_groupSize;
  std::vector<std::string> m_FID, m_SID;

  std::vector<std::string> m_chromosome, m_snpID, m_refAllele, m_altAllele;
  std::vector<int> m_location;
  std::vector<std::vector<double> > m_altFreq, m_maf, m_avgCall, m_rSq;

  CGeneticDataReader *m_geneticDataReader;
  std::vector<double> m_dosage, m_p0, m_p1, m_p2;
  int m_currentSNP;

  int ProcessString(const std::string &_dataString, std::vector<std::string> &_stringsToFill);
  int ReadSNPAdditionalInfo(std::vector<std::vector<double> > &_infotoRead);
  int MissingSNPAdditionalInfo(std::vector<std::vector<double> > &_infotoRead);

  CBDoseReader(const std::string &_filename);
public:
  virtual ~CBDoseReader();

  virtual int GetIndices() = 0;

  bool GetFirst();
  bool GetNext();
  bool GetSNP(int n);

  bool good() const { return m_good; }
  int Format() const { return m_format;  }
  int Version() const { return m_version; }
  int NumSamples() const { return m_SID.size(); }
  int NumSNPs() const { return m_chromosome.size(); }
  int NumGroups() const { return m_groupSize.size(); }
  const std::vector<int> &GroupSize() const { return m_groupSize; }
  const std::vector<std::string> &FamilyID() const { return m_FID; }
  const std::vector<std::string> &SampleID() const { return m_SID; }
  const std::vector<std::string> &Chromosome() const { return m_chromosome; }
  const std::vector<std::string> &SNPID() const { return m_snpID; }
  const std::vector<int> &Location() const { return m_location; }
  const std::vector<std::string> &ReferenceAllele() const { return m_refAllele; }
  const std::vector<std::string> &AlternateAllele() const { return m_altAllele; }
  const std::vector<std::vector<double> > &AlternateAlleleFreq() const { return m_altFreq; }
  const std::vector<std::vector<double> > &MinorAlleleFreq() const { return m_maf; }
  const std::vector<std::vector<double> > &AvgCall() const { return m_avgCall; }
  const std::vector<std::vector<double> > &RSquared() const { return m_rSq; }
  const std::vector<double> &Dosage() const { return m_dosage; }
  const std::vector<double> &P0() const { return m_p0; }
  const std::vector<double> &P1() const { return m_p1; }
  const std::vector<double> &P2() const { return m_p2; }
  const std::vector<int> &Indices() const { return m_snpIndex; }
  int CurrentSNP() const { return m_currentSNP + 1; }
};

class CBDoseReader1 : public CBDoseReader {
protected:
public:
  CBDoseReader1(const std::string &_filename, const std::vector<std::string> &_FID, std::vector<std::string> &_SID,
                const std::vector<std::string> &_SNPID, const std::vector<std::string> &_chromosome,
                const std::vector<int> &_location, const std::vector<std::string> &_refAllele,
                const std::vector<std::string> &_altAllele);
  virtual ~CBDoseReader1() {}

  virtual int GetIndices();
};

class CBDoseReader4 : public CBDoseReader {
protected:
public:
  CBDoseReader4(const std::string &_filename);
  virtual ~CBDoseReader4() {}

  virtual int GetIndices();
};

#endif
