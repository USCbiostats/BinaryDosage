#ifndef BINARYDOSAGEWRITER_H
#define BINARYDOSAGEWRITER_H 1

#ifndef GENETICDATAWRITER_H
#include "GeneticDataWriter.h"
#endif

int WriteFamFile(const std::string &_filename, const std::vector<std::string> &_FID, const std::vector<std::string> &_SID);
int WriteMapFile(const std::string &_filename, const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
	const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele);

class CBinaryDosageWriter {
protected:
	bool m_good;
	std::string m_filename;
	std::fstream m_outfile;
	const int m_version, m_subversion;

	std::streampos m_startDosageData;

	CGeneticDataWriter *m_geneticDataWriter;

	CBinaryDosageWriter(const std::string &_filename, const int _version, const int _subversion);
public:
	virtual ~CBinaryDosageWriter();

	int WriteGeneticData(const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2);
	virtual int UpdateAlternateAlleleFreq(const std::vector<std::vector<double> > &_aaf) { return 0; }

	bool good() const { return m_good; }
	int Version() const { return m_version; }
	int SubVersion() const { return m_subversion; }
};

class CBinaryDosageWriter1 : public CBinaryDosageWriter {
protected:
public:
	CBinaryDosageWriter1(const std::string &_filename, const int _version, const int _subversion,
		const std::string &_famFilename, const std::string &_mapFilename,
		const std::vector<std::string> &_FID, const std::vector<std::string> &_SID,
		const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
		const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
		const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
		const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq) : CBinaryDosageWriter(_filename, _version, _subversion) {}
	virtual ~CBinaryDosageWriter1() {};
};

class CBinaryDosageWriter4 : public CBinaryDosageWriter {
protected:
	std::streampos m_startSamples, m_startSNPs;

	int ProcessString(const std::vector<std::string> &_stringsToWrite, unsigned int &_strLength);
	int ProcessExtraSNPData(const std::vector<std::vector<double> > &_dataToWrite, const unsigned int _groupSize);

	unsigned int GetSNPOptions(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
		const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
		const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
		const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq);
public:
	CBinaryDosageWriter4(const std::string &_filename, const int _version, const int _subversion,
		const std::vector<unsigned int> &_groups, const std::vector<std::string> &_FID, const std::vector<std::string> &_SID,
		const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
		const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
		const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
		const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq);
	virtual ~CBinaryDosageWriter4() {};

	virtual int UpdateAlternateAlleleFreq(const std::vector<std::vector<double> > &_aaf);
};

#endif
