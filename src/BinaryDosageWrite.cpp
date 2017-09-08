#define RProject 1
#ifdef RProject
#include <RcppArmadillo.h>
#define CERR Rcpp::Rcerr
#else
#define CERR std::cerr
#define NA_REAL NAN
#endif
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
#include "BinaryDosageWrite.h"

// ********************************************************************************************
//
// Converting binary dosage file format 4
//
// ********************************************************************************************
CWriteBinaryDosage4::CWriteBinaryDosage4() {
	m_version = 4;
	m_subversion = 2;
	m_numSubjects = 0;
	m_numSNPs = 0;
	m_numGroups = 0;
	m_subjectOptions = 0;
	m_SNPOptions = 0;
	m_subjectStart = 0;
	m_SNPStart = 0;
	m_SNPNameSize = 0;
	m_chromosomeSize = 0;
	m_refAlleleSize = 0;
	m_altAlleleSize = 0;
	m_dosageStart = 0;
	m_subjectIDSize = 0;
	m_familyIDSize = 0;
}

CWriteBinaryDosage4::~CWriteBinaryDosage4() {
	// Close the files - not much else to do
	m_infile.close();
	m_outfile.close();
}
// Write header - stats about the file
int CWriteBinaryDosage4::WriteHeader() {
	char magicWord[8] = { 'b', 'o', 's', 'e', 0x0, 0x0, 0x0, 0x0 };
	switch (m_version) {
	case 1:
		magicWord[5] = 0x1;
		break;
	case 3:
		magicWord[5] = 0x3;
		break;
	case 4:
		magicWord[5] = 0x4;
		break;
	default:
		CERR << "Unknown version number" << std::endl;
		return 1;
	}
	switch (m_subversion) {
	case 1:
		magicWord[7] = 0x1;
		break;
	case 2:
		magicWord[7] = 0x2;
		break;
	default:
		CERR << "Unknown version number" << std::endl;
		return 1;
	}
	m_outfile.write(magicWord, 8);
	m_outfile.write((char *)&m_numSubjects, sizeof(unsigned int));
	m_outfile.write((char *)&m_numSNPs, sizeof(unsigned int));
	m_outfile.write((char *)&m_numGroups, sizeof(unsigned int));
	m_outfile.write((char *)&m_subjectOptions, sizeof(unsigned int));
	m_outfile.write((char *)&m_SNPOptions, sizeof(unsigned int));
	m_subjectStart = 40 + m_numGroups * sizeof(unsigned int);
	//  Rcpp::Rcout << "Subject start:\t" << m_subjectStart << std::endl;
	m_SNPStart = m_subjectStart + m_subjectIDSize + 2 * sizeof(unsigned int);
	//  Rcpp::Rcout << "SNP start:\t" << m_SNPStart << std::endl;
	m_dosageStart = m_SNPStart + 4 * sizeof(unsigned int)+m_SNPNameSize + m_chromosomeSize + m_refAlleleSize + m_altAlleleSize
		+ m_numSNPs * sizeof(unsigned int)+4 * m_numSNPs * sizeof(double);
	//  Rcpp::Rcout << "Dosage start:\t" << m_dosageStart << std::endl;
	m_outfile.write((char *)&m_subjectStart, sizeof(unsigned int));
	m_outfile.write((char *)&m_SNPStart, sizeof(unsigned int));
	m_outfile.write((char *)&m_dosageStart, sizeof(unsigned int));
	m_outfile.write((char *)&m_numSubjects, sizeof(unsigned int));
	return 0;
}
// Write the subject IDs (no family IDs in VCF files)
void CWriteBinaryDosage4::WriteSubjects() {
	unsigned int ui;
	const char tabChar = '\t';

	m_outfile.write((char *)&m_subjectIDSize, sizeof(unsigned int));
	m_outfile.write((char *)&m_familyIDSize, sizeof(unsigned int));
	for (ui = 0; ui < m_numSubjects; ++ui) {
		m_outfile.write(m_subjectID[ui].c_str(), m_subjectID[ui].length());
		m_outfile.write(&tabChar, 1);
	}
}
// Write SNP information - SNP Name, chromosome, location, refAllele, altAllele, altFreq, maf, avgCall, rsq
void CWriteBinaryDosage4::WriteSNPInfo() {
	unsigned int ui;
	const char tabChar = '\t';

	m_outfile.write((char *)&m_SNPNameSize, sizeof(unsigned int));
	m_outfile.write((char *)&m_chromosomeSize, sizeof(unsigned int));
	m_outfile.write((char *)&m_refAlleleSize, sizeof(unsigned int));
	m_outfile.write((char *)&m_altAlleleSize, sizeof(unsigned int));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	for (ui = 0; ui < m_SNPName.size(); ++ui) {
		m_outfile.write(m_SNPName[ui].c_str(), m_SNPName[ui].length());
		m_outfile.write(&tabChar, 1);
	}
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	for (ui = 0; ui < m_chromosome.size(); ++ui) {
		m_outfile.write(m_chromosome[ui].c_str(), m_chromosome[ui].length());
		m_outfile.write(&tabChar, 1);
	}
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	m_outfile.write((char *)&m_location[0], m_location.size() * sizeof(unsigned int));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	for (ui = 0; ui < m_refAllele.size(); ++ui) {
		m_outfile.write(m_refAllele[ui].c_str(), m_refAllele[ui].length());
		m_outfile.write(&tabChar, 1);
	}
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	for (ui = 0; ui < m_altAllele.size(); ++ui) {
		m_outfile.write(m_altAllele[ui].c_str(), m_altAllele[ui].length());
		m_outfile.write(&tabChar, 1);
	}
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	m_outfile.write((char *)&m_altFreq[0], m_altFreq.size() * sizeof(double));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	m_outfile.write((char *)&m_maf[0], m_maf.size() * sizeof(double));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	m_outfile.write((char *)&m_avgCall[0], m_avgCall.size() * sizeof(double));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
	m_outfile.write((char *)&m_rSquared[0], m_rSquared.size() * sizeof(double));
	//  Rcpp::Rcout << m_outfile.tellp() << std::endl;
}

// ********************************************************************************************
//
// Converting VCF to binary dosage file
//
// ********************************************************************************************
// Read in the information file - Little to no error checking
int CWriteBinaryDosageFromVCF::ReadInfoFile() {
	std::string header;
	std::string col4;
	std::string snpN;
	std::string ra;
	std::string aa;
	std::string mja, mna;
	std::string junk;
	bool noMM;
	double af, ma, ac, rs;

	// Just in case something is open
	m_infile.close();
	m_infile.clear();
	// Set SNP options for VCF format
	m_SNPOptions = 0x07f4;
	// Clear out strings
	std::vector<std::string>().swap(m_SNPName);
	std::vector<std::string>().swap(m_chromosome);
	std::vector<unsigned int>().swap(m_location);
	std::vector<std::string>().swap(m_altAllele);
	std::vector<std::string>().swap(m_refAllele);
	std::vector<double>().swap(m_altFreq);
	std::vector<double>().swap(m_maf);
	std::vector<double>().swap(m_avgCall);
	std::vector<double>().swap(m_rSquared);
	m_SNPNameSize = 0;
	m_chromosomeSize = 0;
	m_altAlleleSize = 0;
	m_refAlleleSize = 0;
	// Open the information file
	m_infile.open(m_infoFilename.c_str());
	// Did the file open?
	if (!m_infile.good()) {
		CERR << "Unable to open information file:\t" << m_infoFilename << std::endl;
		return 1;
	}
	// Read in the header line
	m_infile >> junk >> junk >> junk >> col4;
	noMM = !(col4 == "Major");
	getline(m_infile, header);
	// Read in the first line
	if (noMM) {
		m_infile >> snpN >> ra >> aa >> af >> ma;
	}
	else {
		m_infile >> snpN >> ra >> aa >> mja >> mna >> ma;
		af = (mna == aa) ? ma : 1 - ma;
	}
	m_infile >> std::ws;
	// Check if average call(ac), and r-squared(rs) are missing
	if (m_infile.peek() == '-') {
		ac = NA_REAL;
		rs = NA_REAL;
	}
	else {
		m_infile >> ac >> rs;
	}
	// Read the rest of the line - data not used
	getline(m_infile, junk);
	while (!m_infile.fail()) {
		m_SNPName.push_back(snpN);
		m_SNPNameSize += snpN.length();
		m_refAllele.push_back(ra);
		m_refAlleleSize += ra.length();
		m_altAllele.push_back(aa);
		m_altAlleleSize += aa.length();
		m_altFreq.push_back(af);
		m_maf.push_back(ma);
		m_avgCall.push_back(ac);
		m_rSquared.push_back(rs);
		if (noMM) {
			m_infile >> snpN >> ra >> aa >> af >> ma;
		}
		else {
			m_infile >> snpN >> ra >> aa >> mja >> mna >> ma;
			af = (mna == aa) ? ma : 1 - ma;
		}
		m_infile >> std::ws;
		if (m_infile.peek() == '-') {
			ac = NA_REAL;
			rs = NA_REAL;
		}
		else {
			m_infile >> ac >> rs;
		}
		getline(m_infile, junk);
	}
	// Save some usefull information about memory usage
	m_numSNPs = m_SNPName.size();
	m_SNPNameSize += m_numSNPs;
	m_altAlleleSize += m_numSNPs;
	m_refAlleleSize += m_numSNPs;

	//  Rcpp::Rcout << "Number of SNPs\t" << m_numSNPs << std::endl;
	m_infile.close();
	// Convert SNP names to chromosome and location - This is the HRC VCF format
	return GetChromosomesAndLocations();
}

int CWriteBinaryDosageFromVCF::GetChromosomesAndLocations() {
	unsigned int ui;
	std::istringstream iss;
	std::string rchr;
	unsigned rloc;
	bool oneChromosome = true;

	m_chromosome.reserve(m_numSNPs);
	m_location.reserve(m_numSNPs);

	for (ui = 0; ui < m_numSNPs; ++ui) {
		iss.clear();
		std::replace(m_SNPName[ui].begin(), m_SNPName[ui].end(), ':', '\t');
		iss.str(m_SNPName[ui]);
		iss >> rchr >> rloc;
		if (iss.fail()) {
			CERR << "SNP name is not in <chr>:<location> format" << std::endl;
			return 1;
		}
		m_chromosome.push_back(rchr);
		m_chromosomeSize += rchr.length();
		m_location.push_back(rloc);
		if (rchr != m_chromosome[0])
			oneChromosome = false;
	}
	// Check for one chromosome - set option flag and reduce array to one value
	if (oneChromosome) {
		m_SNPOptions |= 0x0008;
		std::vector<std::string>().swap(m_chromosome);
		m_chromosome.push_back(rchr);
		m_chromosomeSize = rchr.length() + 1;
	}
	else {
		m_chromosomeSize += m_numSNPs;
	}
	// Clear out SNPName vector - no longer needed
	std::vector<std::string>().swap(m_SNPName);
	m_SNPNameSize = 0;
	return 0;
}

int CWriteBinaryDosageFromVCF::ReadSubjectInfo() {
	std::string readLine;
	std::string columnName;
	unsigned int ui;
	std::istringstream iss;
	std::string riid;

	// Set subject options for VCF format - no family IDs
	m_subjectOptions = 0;
	// Clear out the subject and family IDs
	m_subjectIDSize = 0;
	m_familyIDSize = 0;
	std::vector<std::string>().swap(m_subjectID);
	std::vector<std::string>().swap(m_familyID);

	// The twelfth line should have the subject IDs
	for (ui = 0; ui < 12; ++ui)
		getline(m_infile, readLine);
	if (readLine[0] != '#' || readLine[1] != 'C') {
		CERR << "Twelfth line of the vcf file does not have the subject IDs" << std::endl;
		return 1;
	}

	// Skip the leading # mark
	iss.str(readLine.substr(1));
	// Read in the column names - currently not used
	for (ui = 0; ui < 9; ++ui)
		iss >> columnName;
	if (iss.fail())
		return 1;

	iss >> riid;
	while (!iss.fail()) {
		m_subjectID.push_back(riid);
		m_subjectIDSize += riid.length();
		iss >> riid;
	}

	m_numSubjects = m_subjectID.size();
	m_subjectIDSize += m_numSubjects;
	//  Rcpp::Rcout << "Number of subjects:\t" << m_numSubjects << std::endl;

	return 0;
}

// Read and write the dosage values
int CWriteBinaryDosageFromVCF::WriteDosages() {
	const short addVal[3][10] = { { 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000 }, { 0, 100, 200, 300, 400, 500, 600, 700, 800, 900 }, { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 } };
	bool oneChr;
	unsigned int numSNPs;
	unsigned int numMissing;
	bool atEOF;
	std::string rChr, rSNPName, rRefAllele, rAltAllele, junk, readValues;
	bool missingSNP;
	unsigned int rLoc;
	unsigned int ui, uj, uk, um;
	unsigned int additionalInfoSize;
	unsigned int dosageSize;
	unsigned int extra1, extra2;
	std::vector<short> sdosage;
	std::vector<short> extraValue;
	short p[4];
	const unsigned int uiZero = 0;
	int pSum, pDose;
	char *aPtr;

	sdosage.resize(m_numSubjects);
	extraValue.resize(3 * m_numSubjects);

	oneChr = (m_chromosome.size() == 1);
	numSNPs = 0;
	numMissing = 0;
	atEOF = false;
	m_infile >> rChr >> rLoc >> rSNPName >> rRefAllele >> rAltAllele >> junk >> junk >> junk >> junk;
	m_infile >> std::ws;
	getline(m_infile, readValues);
	for (ui = 0; ui < m_numSNPs; ++ui) {
		if (!atEOF) {
			if (oneChr)
				missingSNP = !(rChr == m_chromosome[0] && rLoc == m_location[ui] && rRefAllele == m_refAllele[ui] && rAltAllele == m_altAllele[ui]);
			else
				missingSNP = !(rChr == m_chromosome[ui] && rLoc == m_location[ui] && rRefAllele == m_refAllele[ui] && rAltAllele == m_altAllele[ui]);
			if (missingSNP) {
				if (m_subversion == 2) {
					m_outfile.write((char *)&uiZero, sizeof(unsigned int));
				}
				else {
					std::fill(sdosage.begin(), sdosage.end(), 20001);
					m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
				}
				++numMissing;
			}
			else {
				aPtr = &readValues[4];
				additionalInfoSize = 0;
				extra1 = 0;
				extra2 = 0;
				for (uj = 0; uj < m_numSubjects; ++uj, aPtr += 4) {
					memset(p, 0, sizeof(p));
					for (uk = 0; uk < 4; ++uk, ++aPtr) {
						switch (*aPtr) {
						case '0':
							break;
						case '1':
							p[uk] = 10000;
							break;
						case'2':
							p[uk] = 20000;
							break;
						case '.':
							if (uk == 0) {
								p[uk] = 20001;
							}
							else {
								CERR << "VCF format error" << std::endl;
								return 1;
							}
							break;
						default:
							CERR << "VCF format error" << std::endl;
							return 1;
						}
						++aPtr;
						++aPtr;
						if (p[uk] == 20001) {
							++aPtr;
							break;
						}
						for (um = 0; um < 3; ++um, ++aPtr) {
							switch (*aPtr) {
							case '0':
								break;
							case '1':
							case '2':
							case '3':
							case '4':
							case '5':
							case '6':
							case '7':
							case '8':
							case '9':
								p[uk] += addVal[um][*aPtr - '0'];
								break;
							default:
								CERR << "VCF format error" << std::endl;
								return 1;
							}
						}
					}
					sdosage[uj] = p[0];
					if (p[0] == 20001 || m_subversion == 1)
						continue;
					pSum = p[1] + p[2] + p[3];
					pDose = p[2] + p[3] + p[3];
					if (pDose != p[0] || pSum != 10000) {
						sdosage[uj] |= 0x8000;
						extraValue[additionalInfoSize] = p[2] | 0x8000;
						++additionalInfoSize;
						extraValue[additionalInfoSize] = p[1];
						++additionalInfoSize;
						extraValue[additionalInfoSize] = p[3];
						++additionalInfoSize;
						++extra2;
					}
					else if (p[1] != 0 && p[3] != 0) {
						//            if (ui == 2)
						//              Rcpp::Rcout << p[0] << '\t' << p[1] << '\t' << p[2] << '\t' << p[3] << '\t' << pDose << '\t' << pSum << std::endl;
						sdosage[uj] |= 0x8000;
						extraValue[additionalInfoSize] = p[2];
						++additionalInfoSize;
						++extra1;
					}
				}
				//        Rcpp::Rcout << extra1 << '\t' << extra2 << std::endl;
				++numSNPs;
				dosageSize = m_numSubjects + m_numSubjects + additionalInfoSize + additionalInfoSize;
				if (m_subversion == 2)
					m_outfile.write((char *)&dosageSize, sizeof(unsigned int));
				m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
				if (additionalInfoSize > 0)
					m_outfile.write((char *)&extraValue[0], additionalInfoSize + additionalInfoSize);
				m_infile >> rChr;
				if (m_infile.fail()) {
					if (m_infile.eof()) {
						atEOF = true;
					}
					else {
						CERR << "Unknown error when reading VCF file" << std::endl;
						return 1;
					}
				}
				else {
					m_infile >> rLoc >> rSNPName >> rRefAllele >> rAltAllele >> junk >> junk >> junk >> junk;
					m_infile >> std::ws;
					if (m_infile.fail()) {
						CERR << "Error reading SNP information in VCF file" << std::endl;
						return 1;
					}
					else {
						getline(m_infile, readValues);
					}
				}
			}
		}
		else {
			if (m_subversion == 2) {
				m_outfile.write((char *)&uiZero, sizeof(unsigned int));
			}
			else {
				std::fill(sdosage.begin(), sdosage.end(), 20001);
				m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
			}
			++numMissing;
		}
	}
	m_outfile.close();
	//  Rcpp::Rcout << "Number of SNPs\t" << numSNPs << std::endl << "Number missing from VCF file\t" << numMissing << std::endl;
	if (!atEOF) {
		getline(m_infile, readValues);
		if (readValues == "")
			return 0;
		CERR << "Did not reached EOF" << std::endl;
		CERR << "Possible extra SNPs or blank lines at end of VCF file" << std::endl;
		return 1;
	}
	return 0;
}


// Need to change the first two elements to std::vector<string> to allow for mulitple input files
int CWriteBinaryDosageFromVCF::WriteBinaryDosageFile(const std::string &vcfFilename, const std::string &infoFilename, const std::string &outputFilename,
	unsigned int version, unsigned int subversion) {

	m_version = version;
	m_subversion = subversion;
	m_numGroups = 1;
	// Read the information from the info file
	m_infoFilename = infoFilename;
	if (ReadInfoFile())
		return 1;
	// Close and clear m_infile
	m_infile.close();
	m_infile.clear();
	// Open the vcf file
	m_infile.open(vcfFilename.c_str());
	if (!m_infile.good()) {
		CERR << "Unable to open dosage file:\t" << vcfFilename << std::endl;
		return 1;
	}
	// Read in the subject IDs
	if (ReadSubjectInfo())
		return 1;
	// With subject and SNP information read it is now time to write the header
	// Open the output file
	m_outfile.open(outputFilename.c_str(), std::ios::out | std::ios::binary);
	if (!m_outfile.good()) {
		CERR << "Unable to open output file:\t" << outputFilename << std::endl;
		return 1;
	}
	// Write the header - stats about the file
	if (WriteHeader())
		return 1;
	// Write the subject data - family and subject IDs
	WriteSubjects();
	// Write the information about the SNPs
	WriteSNPInfo();
	// Read and write the dosages
	return WriteDosages();
}

// ********************************************************************************************
//
// Converting VCF that has been converted to Impute 2 format to binary dosage file
//
// ********************************************************************************************
int CWriteBinaryDosageFromVCF2Impute::ReadSubjectInfo() {
	std::ifstream infile;
	std::string subName;

	infile.open(m_sampleFilename.c_str());
	// Set subject options for VCF format - no family IDs
	m_subjectOptions = 0;
	// Clear out the subject and family IDs
	m_subjectIDSize = 0;
	m_familyIDSize = 0;
	std::vector<std::string>().swap(m_subjectID);
	std::vector<std::string>().swap(m_familyID);

	infile >> subName;
	while (!infile.fail()) {
		m_subjectID.push_back(subName);
		m_subjectIDSize += subName.length();
		infile >> subName;
	}

	m_numSubjects = m_subjectID.size();
	m_subjectIDSize += m_numSubjects;

	return 0;
}

int CWriteBinaryDosageFromVCF2Impute::WriteDosages() {
	const short addVal[3][10] = { { 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000 }, { 0, 100, 200, 300, 400, 500, 600, 700, 800, 900 }, { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 } };
	bool oneChr;
	unsigned int numSNPs;
	unsigned int numMissing;
	bool atEOF;
	std::string rChr, rSNPName, rRefAllele, rAltAllele, junk, readValues;
	bool missingSNP;
	unsigned int rLoc;
	unsigned int ui, uj, uk, um;
	unsigned int additionalInfoSize;
	unsigned int dosageSize;
	unsigned int extra1, extra2;
	std::vector<short> sdosage;
	std::vector<short> extraValue;
	short p[4];
	const unsigned int uiZero = 0;
	int pSum, pDose;
	char *aPtr;

	sdosage.resize(m_numSubjects);
	extraValue.resize(3 * m_numSubjects);

	oneChr = (m_chromosome.size() == 1);
	numSNPs = 0;
	numMissing = 0;
	atEOF = false;
	m_infile >> rSNPName >> rLoc >> rRefAllele >> rAltAllele;
	m_infile >> std::ws;
	getline(m_infile, readValues);
	for (ui = 0; ui < m_numSNPs; ++ui) {
		if (!atEOF) {
			missingSNP = false;
			if (rLoc != m_location[ui] || rRefAllele != m_refAllele[ui] || rAltAllele != m_altAllele[ui]) {
				missingSNP = true;
			}
			else {
				uj = oneChr ? 0 : ui;
				if (std::strncmp(rSNPName.c_str(), m_chromosome[uj].c_str(), m_chromosome[uj].length()) || rSNPName[m_chromosome[uj].length()] != ':')
					missingSNP = true;
			}
			if (missingSNP) {
				if (m_subversion == 2) {
					m_outfile.write((char *)&uiZero, sizeof(unsigned int));
				}
				else {
					std::fill(sdosage.begin(), sdosage.end(), 20001);
					m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
				}
				++numMissing;
			}
			else {
				aPtr = &readValues[0];
				additionalInfoSize = 0;
				extra1 = 0;
				extra2 = 0;
				for (uj = 0; uj < m_numSubjects; ++uj) {
					memset(p, 0, sizeof(p));
					for (uk = 1; uk < 4; ++uk, ++aPtr) {
						switch (*aPtr) {
						case '0':
							break;
						case '1':
							p[uk] = 10000;
							break;
						case'2':
							p[uk] = 20000;
							break;
						case '.':
							if (uk == 1) {
								p[uk] = 20001;
							}
							else {
								CERR << "VCF format error" << std::endl;
								return 1;
							}
							break;
						default:
							CERR << "VCF format error" << std::endl;
							return 1;
						}
						++aPtr;
						++aPtr;
						if (p[uk] == 20001) {
							aPtr += 3;
							break;
						}
						for (um = 0; um < 3; ++um, ++aPtr) {
							switch (*aPtr) {
							case '0':
								break;
							case '1':
							case '2':
							case '3':
							case '4':
							case '5':
							case '6':
							case '7':
							case '8':
							case '9':
								p[uk] += addVal[um][*aPtr - '0'];
								break;
							default:
								CERR << "VCF format error" << std::endl;
								return 1;
							}
						}
						aPtr += 3;
					}
					if (p[1] == 20001) {
						p[0] = 20001;
					}
					else {
						p[0] = p[2] + p[3] + p[3];
						if (p[0] > 20000)
							p[0] = 20000;
					}
					sdosage[uj] = p[0];
					if (p[0] == 20001 || m_subversion == 1)
						continue;
					pSum = p[1] + p[2] + p[3];
					pDose = p[2] + p[3] + p[3];
					if (pDose != p[0] || pSum != 10000) {
						sdosage[uj] |= 0x8000;
						extraValue[additionalInfoSize] = p[2] | 0x8000;
						++additionalInfoSize;
						extraValue[additionalInfoSize] = p[1];
						++additionalInfoSize;
						extraValue[additionalInfoSize] = p[3];
						++additionalInfoSize;
						++extra2;
					}
					else if (p[1] != 0 && p[3] != 0) {
						//            if (ui == 2)
						//              Rcpp::Rcout << p[0] << '\t' << p[1] << '\t' << p[2] << '\t' << p[3] << '\t' << pDose << '\t' << pSum << std::endl;
						sdosage[uj] |= 0x8000;
						extraValue[additionalInfoSize] = p[2];
						++additionalInfoSize;
						++extra1;
					}
				}
				//        Rcpp::Rcout << extra1 << '\t' << extra2 << std::endl;
				++numSNPs;
				dosageSize = m_numSubjects + m_numSubjects + additionalInfoSize + additionalInfoSize;
				if (m_subversion == 2)
					m_outfile.write((char *)&dosageSize, sizeof(unsigned int));
				m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
				if (additionalInfoSize > 0)
					m_outfile.write((char *)&extraValue[0], additionalInfoSize + additionalInfoSize);
				m_infile >> rSNPName;
				if (m_infile.fail()) {
					if (m_infile.eof()) {
						atEOF = true;
					}
					else {
						CERR << "Unknown error when reading VCF file" << std::endl;
						return 1;
					}
				}
				else {
					m_infile >> rLoc >> rRefAllele >> rAltAllele;
					m_infile >> std::ws;
					if (m_infile.fail()) {
						CERR << "Error reading SNP information in VCF file" << std::endl;
						return 1;
					}
					else {
						getline(m_infile, readValues);
					}
				}
			}
		}
		else {
			if (m_subversion == 2) {
				m_outfile.write((char *)&uiZero, sizeof(unsigned int));
			}
			else {
				std::fill(sdosage.begin(), sdosage.end(), 20001);
				m_outfile.write((char *)&sdosage[0], m_numSubjects + m_numSubjects);
			}
			++numMissing;
		}
	}
	m_outfile.close();
	//  Rcpp::Rcout << "Number of SNPs\t" << numSNPs << std::endl << "Number missing from VCF file\t" << numMissing << std::endl;
	if (!atEOF) {
		getline(m_infile, readValues);
		if (readValues == "")
			return 0;
		CERR << "Did not reached EOF" << std::endl;
		CERR << "Possible extra SNPs or blank lines at end of VCF file" << std::endl;
		return 1;
	}
	return 0;
}

int CWriteBinaryDosageFromVCF2Impute::WriteBinaryDosageFile(const std::string &vcfFilename, const std::string &infoFilename, const std::string &sampleFilename,
	const std::string &outputFilename, unsigned int version, unsigned int subversion) {
	m_sampleFilename = sampleFilename;
	return CWriteBinaryDosageFromVCF::WriteBinaryDosageFile(vcfFilename, infoFilename, outputFilename, version, subversion);
}
