#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "GeneticDataWriter.h"
#include "BinaryDosageWriter.h"

enum class Header4pos {
	header = 0,
	version = 4,
	numSub = 8,
	numSNPs = 12,
	numGroups = 16,
	subOptions = 20,
	snpOptions = 24,
	startSub = 28,
	startSNP = 32,
	startDosage = 36,
	startGroups = 40
};

int WriteFamFile(const std::string &_filename, const std::vector<std::string> &_FID, const std::vector<std::string> &_SID) {
	std::ofstream outfile;
	std::vector<std::string>::const_iterator strIt1, strIt2;

	if (_SID.size() == 0)
		return 1;
	if (_FID.size() != 0 && _FID.size() != _SID.size())
		return 1;

	outfile.open(_filename.c_str());
	if (!outfile.good())
		return 1;

	if (_FID.size() == 0) {
		for (strIt1 = _SID.begin(); strIt1 != _SID.end(); ++strIt1)
			outfile << *strIt1 << "\t0\t0\t9\t9" << std::endl;
	}
	else {
		strIt2 = _FID.begin();
		for (strIt1 = _SID.begin(); strIt1 != _SID.end(); ++strIt1, ++strIt2)
			outfile << *strIt1 << '\t' << *strIt2 << "\t0\t0\t9\t9" << std::endl;
	}
	outfile.close();
	return 0;
}

int WriteMapFile(const std::string &_filename, const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
	const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele) {
	unsigned int numSNPs;
	std::ofstream outfile;
	std::vector<std::string>::const_iterator strIt1, strIt2, strIt3, strIt4;
	std::vector<int>::const_iterator dblIt;

	numSNPs = _chromosome.size();
	if (numSNPs == 0)
		return 1;
	if (_snpID.size() != 0 && _snpID.size() != numSNPs)
		return 1;
	if (_location.size() != numSNPs || _refAllele.size() != numSNPs || _altAllele.size() != numSNPs)
		return 1;

	outfile.open(_filename.c_str());
	if (!outfile.good())
		return 1;

	if (_snpID.size() == 0) {
		strIt3 = _refAllele.begin();
		strIt4 = _refAllele.begin();
		dblIt = _location.begin();
		for (strIt1 = _chromosome.begin(); strIt1 != _chromosome.end(); ++strIt1, ++strIt3, ++strIt4, ++dblIt)
			outfile << *strIt1 << '\t' << *strIt1 << ':' << *dblIt << "\t0\t" << *dblIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
	}
	else {
		strIt2 = _snpID.begin();
		strIt3 = _refAllele.begin();
		strIt4 = _refAllele.begin();
		dblIt = _location.begin();
		for (strIt1 = _chromosome.begin(); strIt1 != _chromosome.end(); ++strIt1, ++strIt2, ++strIt3, ++strIt4, ++dblIt)
			outfile << *strIt1 << '\t' << *strIt2 << "\t0\t" << *dblIt << '\t' << *strIt3 << '\t' << *strIt4 << std::endl;
	}
	outfile.close();
	return 0;
}

CBinaryDosageWriter::CBinaryDosageWriter(const std::string &_filename, const int _version, const int _subversion) : m_version(_version), m_subversion(_subversion) {
  std::ofstream outfile;

	m_filename = _filename;
	m_geneticDataWriter = NULL;

	outfile.open(m_filename.c_str());
	if (!outfile.good()) {
	  outfile.close();
	  m_good = false;
	  return;
	}
	outfile.close();

	m_outfile.open(m_filename.c_str(), std::ios_base::out | std::ios_base::in | std::ios_base::binary | std::ios_base::trunc);
	if (!m_outfile.good()) {
		m_good = false;
		return;
	}
	m_good = true;
}

CBinaryDosageWriter::~CBinaryDosageWriter() {
	if (m_geneticDataWriter)
		delete m_geneticDataWriter;
	m_outfile.close();
}

int CBinaryDosageWriter::WriteGeneticData(const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
	if (!m_good)
		return 1;
	if (m_geneticDataWriter->WriteData(m_outfile, _dosage, _p0, _p1, _p2)) {
		m_good = false;
		return 1;
	}
	return 0;
}
CBinaryDosageWriter4::CBinaryDosageWriter4(const std::string &_filename, const int _version, const int _subversion,
	const std::vector<unsigned int> &_groups, const std::vector<std::string> &_FID, const std::vector<std::string> &_SID,
	const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
	const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
	const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
	const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq) : CBinaryDosageWriter(_filename, _version, _subversion) {
	unsigned int numSamples, numSNPs, numGroups;
	int startSamples, startSNPs, startDosages;
	unsigned int sampleOptions, snpOptions;
	unsigned int sidLength, fidLength;
	unsigned int snpLength, chrLength, refLength, altLength;
	const char dosageHeader[4] = { 'b', 'o', 's', 'e' };
	const char zeroChar = 0;
	const int zeroInt = 0;
	char cver, csubver;
	std::string s1;
	std::vector<unsigned int>::const_iterator usIt;

	if (!m_good)
		return;
	m_good = false;

	numGroups = _groups.size();
	if (numGroups == 0)
		return;

	numSamples = 0;
	for (usIt = _groups.begin(); usIt != _groups.end(); ++usIt)
		numSamples += *usIt;
	if (numSamples == 0)
		return;
	if (_SID.size() != numSamples)
		return;
	if (_FID.size() != 0 && _FID.size() != numSamples)
		return;

	numSNPs = _chromosome.size();
	if (numSNPs == 0)
		return;
	if (_snpID.size() != 0 && _snpID.size() != numSNPs)
		return;
	if (_location.size() != numSNPs || _refAllele.size() != numSNPs || _altAllele.size() != numSNPs)
		return;
	if (_altFreq.size() != 0 && _altFreq.size() != numSNPs)
		return;
	if (_maf.size() != 0 && _maf.size() != numSNPs)
		return;
	if (_avgCall.size() != 0 && _avgCall.size() != numSNPs)
		return;
	if (_rSq.size() != 0 && _rSq.size() != numSNPs)
		return;

	cver = (char)m_version;
	csubver = (char)m_subversion;
	m_outfile.write(dosageHeader, 4);
	m_outfile.write(&zeroChar, 1);
	m_outfile.write(&cver, 1);
	m_outfile.write(&zeroChar, 1);
	m_outfile.write(&csubver, 1);

	m_outfile.write((char *)&numSamples, sizeof(unsigned int));
	m_outfile.write((char *)&numSNPs, sizeof(unsigned int));
	m_outfile.write((char *)&numGroups, sizeof(unsigned int));

	if (_FID.size() == 0 || _FID[0] == "")
		sampleOptions = 1;
	else
		sampleOptions = 0;
	m_outfile.write((char *)&sampleOptions, sizeof(unsigned int));
	snpOptions = GetSNPOptions(_chromosome, _snpID, _location, _refAllele, _altAllele, _altFreq, _maf, _avgCall, _rSq);
	m_outfile.write((char *)&snpOptions, sizeof(unsigned int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)_groups.data(), numGroups * sizeof(unsigned int));

	m_startSamples = m_outfile.tellp();
	startSamples = (int)m_startSamples;
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	ProcessString(_SID, sidLength);
	if (_FID.size() == 0 || _FID[0] == "")
		fidLength = 0;
	else
		ProcessString(_FID, fidLength);

	m_startSNPs = m_outfile.tellp();
	startSNPs = (int)m_startSNPs;
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	m_outfile.write((char *)&zeroInt, sizeof(int));
	if (_snpID.size() == 0 || _snpID[0] == "")
	  snpLength = 0;
	else
	  ProcessString(_snpID, snpLength);
	if (snpOptions & 0x0008) {
		s1 = _chromosome[0];
		m_outfile.write(s1.c_str(), s1.length());
		m_outfile.write(&zeroChar, 1);
		chrLength = s1.length() + 1;
	}
	else {
		ProcessString(_chromosome, chrLength);
	}
	m_outfile.write((char *)_location.data(), numSNPs * sizeof(unsigned int));
	if (_refAllele.size() != 0)
		ProcessString(_refAllele, refLength);
	else
		refLength = 0;
	if (_altAllele.size() != 0)
		ProcessString(_altAllele, altLength);
	else
		altLength = 0;

	if (ProcessExtraSNPData(_altFreq, numGroups))
		return;
	if (ProcessExtraSNPData(_maf, numGroups))
		return;
	if (ProcessExtraSNPData(_avgCall, numGroups))
		return;
	if (ProcessExtraSNPData(_rSq, numGroups))
		return;
	m_startDosageData = m_outfile.tellp();
	startDosages = (int)m_startDosageData;

	m_outfile.seekp((int)Header4pos::startSub);
	m_outfile.write((char *)&startSamples, sizeof(int));
	m_outfile.write((char *)&startSNPs, sizeof(int));
	m_outfile.write((char *)&startDosages, sizeof(int));

	m_outfile.seekp(m_startSamples);
	m_outfile.write((char *)&sidLength, sizeof(int));
	m_outfile.write((char *)&fidLength, sizeof(int));

	m_outfile.seekp(m_startSNPs);
	m_outfile.write((char *)&snpLength, sizeof(int));
	m_outfile.write((char *)&chrLength, sizeof(int));
	m_outfile.write((char *)&refLength, sizeof(int));
	m_outfile.write((char *)&altLength, sizeof(int));

	m_outfile.seekp(m_startDosageData);
	if (m_outfile.good())
		m_good = true;

	if (m_subversion == 1)
		m_geneticDataWriter = new CDosageDataWriter(10000, numSamples);
	else
		m_geneticDataWriter = new CGeneticDataWriter3(10000, numSamples);
}

unsigned int CBinaryDosageWriter4::GetSNPOptions(const std::vector<std::string> &_chromosome, const std::vector<std::string> &_snpID,
	const std::vector<int> &_location, const std::vector<std::string> &_refAllele, const std::vector<std::string> &_altAllele,
	const std::vector<std::vector<double> > &_altFreq, const std::vector<std::vector<double> > &_maf,
	const std::vector<std::vector<double> > &_avgCall, const std::vector<std::vector<double> > &_rSq) {
	unsigned int retVal = 0x0014;
	std::string s1;
	std::vector<std::string>::const_iterator strIt;

	if (_snpID.size() != 0)
		retVal |= 0x0002;
	s1 = _chromosome[0];
	for (strIt = _chromosome.begin(); strIt != _chromosome.end(); ++strIt) {
		if (*strIt != s1)
			break;
	}
	if (strIt == _chromosome.end())
		retVal |= 0x0008;
	if (_refAllele.size() != 0)
		retVal |= 0x0020;
	if (_altAllele.size() != 0)
		retVal |= 0x0040;
	if (_altFreq.size() != 0)
		retVal |= 0x0080;
	if (_maf.size() != 0)
		retVal |= 0x0100;
	if (_avgCall.size() != 0)
		retVal |= 0x0200;
	if (_rSq.size() != 0)
		retVal |= 0x0400;
	return retVal;
}

int CBinaryDosageWriter4::ProcessString(const std::vector<std::string> &_stringsToWrite, unsigned int &_strLength) {
	std::vector<std::string>::const_iterator strIt;
	std::ostringstream oss;
	std::string s1;

	oss.str("");
	for (strIt = _stringsToWrite.begin(); strIt != _stringsToWrite.end(); ++strIt)
		oss << *strIt << '\t';
	s1 = oss.str();
	s1[s1.length() - 1] = 0;
	m_outfile.write(s1.c_str(), s1.length());
	_strLength = s1.length();

	return 0;
}

int CBinaryDosageWriter4::ProcessExtraSNPData(const std::vector<std::vector<double> > &_dataToWrite, const unsigned int _groupSize) {
	std::vector<std::vector<double> >::const_iterator vdblIt;

	for (vdblIt = _dataToWrite.begin(); vdblIt != _dataToWrite.end(); ++vdblIt) {
		if (vdblIt->size() != _groupSize)
			return 1;
		m_outfile.write((char *)vdblIt->data(), _groupSize * sizeof(double));
	}
	return 0;
}

int CBinaryDosageWriter4::UpdateAlternateAlleleFreq(const std::vector<std::vector<double> > &_aaf) {
  int numSNPs, numGroups, snpOptions, chrSize, snpSize, refSize, altSize, jumpSize;
  if (!m_good)
    return 1;

  m_outfile.seekg((int)Header4pos::numSNPs);
  m_outfile.read((char *)&numSNPs, sizeof(int));
  m_outfile.read((char *)&numGroups, sizeof(int));
  if (numSNPs == 0)
    return 1;

  m_outfile.seekg((int)Header4pos::snpOptions);
  m_outfile.read((char *)&snpOptions, sizeof(int));
  if ((snpOptions & 0x0080) == 0)
    return 1;

  m_outfile.seekg(m_startSNPs);
  m_outfile.read((char *)&chrSize, sizeof(int));
  m_outfile.read((char *)&snpSize, sizeof(int));
  m_outfile.read((char *)&altSize, sizeof(int));
  m_outfile.read((char *)&refSize, sizeof(int));

  jumpSize = chrSize + snpSize + altSize + refSize + numSNPs * sizeof(int);
  m_outfile.seekg(jumpSize, std::ios_base::cur);
  m_outfile.seekp(m_outfile.tellg());

  return ProcessExtraSNPData(_aaf, numGroups);
}
