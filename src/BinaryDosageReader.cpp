#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include "GeneticDataReader.h"
#include "BinaryDosageReader.h"

//#define BINARYDOSAGEREADER_DEBUG 1

const char BDHeader[4] = { 'b', 'o', 's', 'e' };
const int BDHeaderSize = 8;

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

int GetBinaryDosageFormat(const std::string &_filename, int &_version, int &_subversion) {
	char header[BDHeaderSize];
	std::ifstream infile;
	int fullVersion;
	int retVal;

	infile.open(_filename.c_str(), std::ios_base::in | std::ios_base::binary);
	if (!infile.good())
		return 1;
	infile.read(header, BDHeaderSize);
	infile.close();
	if (!infile.good())
		return 1;
	if (std::memcmp(header, BDHeader, 4))
		return 1;
	_version = header[5];
	_subversion = header[7];

	fullVersion = _version * 100 + _subversion;
	switch (fullVersion) {
	case 101:
	case 102:
	case 201:
	case 202:
	case 301:
	case 302:
	case 401:
	case 402:
		retVal = 0;
		break;
	default:
		retVal = 1;
		break;
	}
	return retVal;
}

int ReadFamFile(const std::string &_filename, std::vector<std::string> &_FID, std::vector<std::string> &_SID) {
	std::ifstream infile;
	std::istringstream iss;
	std::string s1, s2;
	unsigned int numCol, numSamples;
	unsigned int ui;
	std::vector<std::string>::iterator strIt1, strIt2;
	bool hasFID;

	infile.open(_filename.c_str());
	if (!infile.good())
		return 1;

	std::getline(infile, s1);
	if (!infile.good())
		return 1;
	iss.str(s1);
	numCol = 0;
	do {
		iss >> s2;
		if (iss.fail())
			break;
		++numCol;
	} while (1);
	switch (numCol) {
	case 1:
	case 5:
		hasFID = false;
		break;
	case 2:
	case 6:
		hasFID = true;
		break;
	default:
		infile.close();
		return 1;
		break;
	}

	infile.seekg(0);
	numSamples = 0;
	do {
		std::getline(infile, s1);
		if (infile.eof())
			break;
		++numSamples;
	} while (1);
	_SID.resize(numSamples);
	_FID.resize(numSamples);
	infile.clear();
	infile.seekg(0);
	strIt2 = _FID.begin();
	for (strIt1 = _SID.begin(); strIt1 != _SID.end(); ++strIt1, ++strIt2) {
		std::getline(infile, s1);
		if (!infile.good())
			return 1;
		iss.clear();
		iss.str(s1);
		ui = 1;
		if (hasFID) {
			iss >> s2;
			*strIt2 = s2;
		}
		else {
			*strIt2 = "";
		}
		iss >> s2;
		*strIt1 = s2;
		for (; ui < numCol; ++ui)
			iss >> s2;
		if (iss.fail()) {
			infile.close();
			return 1;
		}
	}

	infile.close();
	return 0;
}

int ReadMapFile(const std::string &_filename, std::vector<std::string> &_chromosome, std::vector<std::string> &_snpID,
	std::vector<int> &_location, std::vector<std::string> &_refAllele, std::vector<std::string> &_altAllele) {
	std::ifstream infile;
	std::istringstream iss;
	std::string s1, s2, s3, s4, s5;
	int x1, x2;
	unsigned int numCol, numSNPs;
	std::vector<std::string>::iterator chromosomeIt, snpIt, refIt, altIt;
	std::vector<int>::iterator locIt;

	infile.open(_filename.c_str());
	if (!infile.good())
		return 1;

	std::getline(infile, s1);
	if (!infile.good())
		return 1;

	iss.str(s1);
	numCol = 0;
	do {
		iss >> s2;
		if (iss.fail())
			break;
		++numCol;
	} while (1);
	if (numCol != 6) {
		infile.close();
		return 1;
	}

	infile.seekg(0);
	numSNPs = 0;
	do {
		std::getline(infile, s1);
		if (infile.eof())
			break;
		++numSNPs;
	} while (1);

	_chromosome.resize(numSNPs);
	_snpID.resize(numSNPs);
	_location.resize(numSNPs);
	_refAllele.resize(numSNPs);
	_altAllele.resize(numSNPs);

	infile.clear();
	infile.seekg(0);
	snpIt = _snpID.begin();
	locIt = _location.begin();
	refIt = _refAllele.begin();
	altIt = _altAllele.begin();
	for (chromosomeIt = _chromosome.begin(); chromosomeIt != _chromosome.end(); ++chromosomeIt, ++snpIt, ++locIt, ++refIt, ++altIt) {
		std::getline(infile, s1);
		if (!infile.good())
			return 1;
		iss.clear();
		iss.str(s1);
		iss >> s2 >> s3 >> x1 >> x2 >> s4 >> s5;
		if (iss.fail()) {
			infile.close();
			return 1;
		}
		*chromosomeIt = s1;
		*snpIt = s2;
		*locIt = x2;
		*refIt = s4;
		*altIt = s5;
	}

	infile.close();
	return 0;
}

CBinaryDosageReader::CBinaryDosageReader(const std::string &_filename) {
	m_good = false;
	m_version = 0;
	m_subversion = 0;
	m_filename = _filename;
	m_geneticDataReader = NULL;
	m_currentSNP = 0;
	if (GetBinaryDosageFormat(m_filename, m_version, m_subversion))
		return;
	m_infile.open(_filename.c_str(), std::ios_base::in | std::ios_base::binary);
	if (!m_infile.good())
		return;
	m_good = true;
}

CBinaryDosageReader::~CBinaryDosageReader() {
	if (m_geneticDataReader)
		delete m_geneticDataReader;
	m_infile.close();
}

int CBinaryDosageReader::ProcessString(const std::string &_dataString, std::vector<std::string> &_stringsToFill) {
	std::vector<std::string>::iterator strIt;
	std::istringstream iss;
	std::string s1;

	iss.str(_dataString);
	for (strIt = _stringsToFill.begin(); strIt != _stringsToFill.end(); ++strIt) {
		iss >> s1;
		*strIt = s1;
	}
	if (iss.fail()) {
		m_good = false;
		return 1;
	}
	return 0;
}

int CBinaryDosageReader::ReadSNPAdditionalInfo(std::vector<std::vector<double> > &_infotoRead) {
	std::vector<std::vector<double> >::iterator vdblIt;

	_infotoRead.resize(NumSNPs());
	for (vdblIt = _infotoRead.begin(); vdblIt != _infotoRead.end(); ++vdblIt)
		vdblIt->resize(NumGroups());
	for (vdblIt = _infotoRead.begin(); vdblIt != _infotoRead.end(); ++vdblIt) {
		m_infile.read((char *)vdblIt->data(), NumGroups() * sizeof(double));
	}
	if (!m_infile.good()) {
		m_good = false;
		return 1;
	}
	return 0;
}

int CBinaryDosageReader::MissingSNPAdditionalInfo(std::vector<std::vector<double> > &_infotoRead) {
	std::vector<std::vector<double> >::iterator vdblIt;
	std::vector<double>::iterator dblIt;

	_infotoRead.resize(NumSNPs());
	for (vdblIt = _infotoRead.begin(); vdblIt != _infotoRead.end(); ++vdblIt)
		vdblIt->resize(NumGroups());
	for (vdblIt = _infotoRead.begin(); vdblIt != _infotoRead.end(); ++vdblIt) {
		for (dblIt = vdblIt->begin(); dblIt != vdblIt->end(); ++dblIt)
			*dblIt = std::numeric_limits<double>::quiet_NaN();
	}
	return 0;
}

bool CBinaryDosageReader::GetFirst() {
	if (!m_good)
		return false;
	m_currentSNP = 0;
	m_infile.seekg(m_startDosageData);
	if (m_geneticDataReader->ReadData(m_infile, m_dosage, m_p0, m_p1, m_p2))
		m_good = false;
	return m_good;
}

bool CBinaryDosageReader::GetNext() {
	if (!m_good)
		return false;
	++m_currentSNP;
	if (m_currentSNP == NumSNPs()) {
		--m_currentSNP;
		return false;
	}
	if (m_geneticDataReader->ReadData(m_infile, m_dosage, m_p0, m_p1, m_p2))
		m_good = false;
	return m_good;
}

bool CBinaryDosageReader::GetSNP(unsigned int n) {
	if (!m_good)
		return false;
	if (n == 0 || n > NumSNPs())
		return false;
	--n;
	if (n == m_currentSNP)
		return true;
	if (n == 0)
		return GetFirst();
	if (n < m_currentSNP) {
		m_infile.seekg(m_startDosageData);
		m_geneticDataReader->SkipSNP(m_infile);
		m_currentSNP = 0;
	}
	while (m_currentSNP < n - 1) {
		m_geneticDataReader->SkipSNP(m_infile);
		++m_currentSNP;
	}
	return GetNext();
}

CBinaryDosageReader4::CBinaryDosageReader4(const std::string &_filename) : CBinaryDosageReader(_filename) {
	int numSamples, numSNPs, numGroups, sampleOptions, snpOptions;
	int startSampleData, startSNPData, startDosageData;
	int sidSize, fidSize;
	int chromosomeSize, snpSize, refSize, altSize;
	int maxSize;
	char *dataString;
	std::string s1;
	std::vector<std::string>::iterator strIt, strIt2;
	std::vector<int>::iterator intIt;

	if (!m_good)
		return;
	m_good = false;

	m_infile.seekg((int)Header4pos::numSub);
	m_infile.read((char *)&numSamples, sizeof(int));
	m_infile.read((char *)&numSNPs, sizeof(int));
	m_infile.read((char *)&numGroups, sizeof(int));
	m_infile.read((char *)&sampleOptions, sizeof(int));
	m_infile.read((char *)&snpOptions, sizeof(int));
	m_infile.read((char *)&startSampleData, sizeof(int));
	m_infile.read((char *)&startSNPData, sizeof(int));
	m_infile.read((char *)&startDosageData, sizeof(int));
#ifdef BINARYDOSAGEREADER_DEBUG
	std::cout << numSamples << '\t' << numSNPs << '\t' << numGroups << '\t'
           << std::hex << sampleOptions << '\t' << snpOptions << '\t'
           << std::dec << startSampleData << '\t' << startSNPData << '\t' << startDosageData << std::endl;
# endif
	m_startDosageData = startDosageData;

	if (!m_infile.good())
		return;

	m_groupSize.resize(numGroups);
	m_infile.read((char *)m_groupSize.data(), numGroups * sizeof(unsigned int));

	if (!m_infile.good())
		return;

	m_infile.read((char *)&sidSize, sizeof(int));
	maxSize = sidSize;
	m_infile.read((char *)&fidSize, sizeof(int));
	if (fidSize > maxSize)
		maxSize = fidSize;
	m_infile.seekg(sidSize + fidSize, std::ios_base::cur);
	m_infile.read((char *)&snpSize, sizeof(int));
	if (snpSize > maxSize)
		maxSize = snpSize;
	m_infile.read((char *)&chromosomeSize, sizeof(int));
	if (chromosomeSize > maxSize)
		maxSize = chromosomeSize;
	m_infile.read((char *)&refSize, sizeof(int));
	if (refSize > maxSize)
		maxSize = refSize;
	m_infile.read((char *)&altSize, sizeof(int));
	if (altSize > maxSize)
		maxSize = altSize;
	if (!m_infile.good())
		return;
#ifdef BINARYDOSAGEREADER_DEBUG
	std::cout << sidSize << '\t' << fidSize << '\t' << snpSize << '\t' << chromosomeSize << '\t'
           << refSize << '\t' << altSize << std::endl;
#endif

	dataString = new char[maxSize + 1];
	m_infile.seekg(startSampleData + 2*sizeof(int));

	m_SID.resize(numSamples);
	m_infile.read(dataString, sidSize);
	dataString[sidSize] = 0;
	s1 = dataString;
	if (ProcessString(s1, m_SID)) {
		if (dataString)
			delete[] dataString;
		return;
	}

	m_FID.resize(numSamples);
	if (fidSize > 0) {
		m_infile.read(dataString, fidSize);
		dataString[fidSize] = 0;
		s1 = dataString;
		if (ProcessString(s1, m_FID)) {
			if (dataString)
				delete[] dataString;
			return;
		}
	}
	else {
		for (strIt = m_FID.begin(); strIt != m_FID.end(); ++strIt)
			*strIt = "";
	}

	m_infile.seekg(startSNPData + 4 * sizeof(int));
	m_snpID.resize(numSNPs);
	if (snpSize > 0) {
		m_infile.read(dataString, snpSize);
		dataString[snpSize] = 0;
		s1 = dataString;
		if (ProcessString(s1, m_snpID)) {
			if (dataString)
				delete[] dataString;
			return;
		}
	}

	m_chromosome.resize(numSNPs);
	m_infile.read(dataString, chromosomeSize);
	dataString[chromosomeSize - 1] = 0;
	s1 = dataString;
	if (snpOptions & 0x0008) {
		for (strIt = m_chromosome.begin(); strIt != m_chromosome.end(); ++strIt)
			*strIt = s1;
	}
	else {
		if (ProcessString(s1, m_chromosome)) {
			if (dataString)
				delete[] dataString;
			return;
		}
	}

	m_location.resize(numSNPs);
	m_infile.read((char *)m_location.data(), numSNPs * sizeof(int));

	if (snpSize == 0) {
	  intIt = m_location.begin();
	  strIt2 = m_chromosome.begin();
	  for (strIt = m_snpID.begin(); strIt != m_snpID.end(); ++strIt, ++strIt2, ++intIt)
	    *strIt = *strIt2 + ":" + std::to_string(*intIt);
	}

	m_refAllele.resize(numSNPs);
	if (refSize > 0) {
		m_infile.read(dataString, refSize);
		dataString[refSize] = 0;
		s1 = dataString;
		if (ProcessString(s1, m_refAllele)) {
			if (dataString)
				delete[] dataString;
			return;
		}
	}
	else {
		for (strIt = m_refAllele.begin(); strIt != m_refAllele.end(); ++strIt)
			*strIt = "1";
	}

	m_altAllele.resize(numSNPs);
	if (altSize > 0) {
		m_infile.read(dataString, altSize);
		dataString[altSize] = 0;
		s1 = dataString;
		if (ProcessString(s1, m_altAllele)) {
			if (dataString)
				delete[] dataString;
			return;
		}
	}
	else {
		for (strIt = m_altAllele.begin(); strIt != m_altAllele.end(); ++strIt)
			*strIt = "2";
	}

	if (dataString)
		delete[] dataString;
	dataString = NULL;

	if (snpOptions & 0x0080) {
		if (ReadSNPAdditionalInfo(m_altFreq))
			return;
	}
	else {
		MissingSNPAdditionalInfo(m_altFreq);
	}

	if (snpOptions & 0x0100) {
		if (ReadSNPAdditionalInfo(m_maf))
			return;
	}
	else {
		MissingSNPAdditionalInfo(m_maf);
	}

	if (snpOptions & 0x0200) {
		if (ReadSNPAdditionalInfo(m_avgCall))
			return;
	}
	else {
		MissingSNPAdditionalInfo(m_avgCall);
	}

	if (snpOptions & 0x0400) {
		if (ReadSNPAdditionalInfo(m_rSq))
			return;
	}
	else {
		MissingSNPAdditionalInfo(m_rSq);
	}


	m_dosage.resize(NumSamples());
	m_p0.resize(NumSamples());
	m_p1.resize(NumSamples());
	m_p2.resize(NumSamples());

	if (m_subversion == 1)
		m_geneticDataReader = new CGeneticDataReader1(10000, NumSamples());
	else
		m_geneticDataReader = new CGeneticDataReader3(10000, NumSamples());

	m_good = true;
	m_good = GetFirst();
}
