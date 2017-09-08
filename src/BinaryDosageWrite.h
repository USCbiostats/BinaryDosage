#ifndef BINARYDOSAGEWRITE_H
#define BINARYDOSAGEWRITE_H 1
// ******************************************************************
// Creating Binary Dosage Files
// ******************************************************************
class CWriteBinaryDosage4 {
protected:
	unsigned int m_version;
	unsigned int m_subversion;
	unsigned int m_numSubjects;
	unsigned int m_numSNPs;
	unsigned int m_numGroups;
	unsigned int m_subjectOptions;
	unsigned int m_SNPOptions;
	unsigned int m_subjectStart;
	unsigned int m_SNPStart;
	unsigned int m_dosageStart;

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

	std::vector<std::string> m_subjectID;
	std::vector<std::string> m_familyID;
	unsigned int m_subjectIDSize;
	unsigned int m_familyIDSize;

	std::ifstream m_infile;
	std::ofstream m_outfile;

	virtual int ReadInfoFile() = 0;
	virtual int ReadSubjectInfo() = 0;
	virtual int WriteHeader();
	virtual void WriteSubjects();
	virtual void WriteSNPInfo();
	virtual int WriteDosages() = 0;
public:
	CWriteBinaryDosage4();
	virtual ~CWriteBinaryDosage4();
};

class CWriteBinaryDosageFromVCF : public CWriteBinaryDosage4 {
protected:
	std::string m_infoFilename;

	virtual int ReadInfoFile();
	int GetChromosomesAndLocations();
	virtual int ReadSubjectInfo();
	virtual int WriteDosages();
public:
	CWriteBinaryDosageFromVCF() : CWriteBinaryDosage4() { m_infoFilename = ""; }
	virtual ~CWriteBinaryDosageFromVCF() {}

	int WriteBinaryDosageFile(const std::string &vcfFilename, const std::string &infoFilename, const std::string &outputFilename,
		unsigned int version = 4, unsigned int subversion = 2);
};

class CWriteBinaryDosageFromVCF2Impute : public CWriteBinaryDosageFromVCF {
protected:
	std::string m_sampleFilename;

	virtual int ReadSubjectInfo();
	virtual int WriteDosages();
public:
	CWriteBinaryDosageFromVCF2Impute() : CWriteBinaryDosageFromVCF() { m_sampleFilename = ""; }
	virtual ~CWriteBinaryDosageFromVCF2Impute() {}

	int WriteBinaryDosageFile(const std::string &vcfFilename, const std::string &infoFilename, const std::string &sampleFilename,
		const std::string &outputFilename, unsigned int version = 4, unsigned int subversion = 2);
};
#endif
