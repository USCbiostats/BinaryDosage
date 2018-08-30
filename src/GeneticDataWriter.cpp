#include <fstream>
#include <cmath>
#include <vector>
#include "GeneticDataWriter.h"

// Tolerance for calculated values for error checking
// Id est, |1. - p0 - p1 - p2| and
// |dosage - p1 - 2p2| must be less than the tolerance
// If not, an error is returned
const double GeneticTolerance = 0.01;
// Tolerance for calculated values for writing values
// Id est, |1. - p0 - p1 - p2| and
// |dosage - p1 - 2p2| must be less than the tolerance
// If not, all four values are written to file
const double WritingTolerance = 0.000001;

////////////////////////////////////////////////////////////////////////////////
//							CGeneticDataWriter
////////////////////////////////////////////////////////////////////////////////
// Base class for writing genetic data to binary dosage files
// When it is created the scale parameter and number of observations
// must be specified. These are constants and cannot be changed.

// Constructor
CGeneticDataWriter::CGeneticDataWriter(const unsigned short _scale, const unsigned int _sampleSize) : m_scale(_scale), m_sampleSize(_sampleSize) {
	// The scale is changed to a double. This helps the speed of
	// execution by not having to constantly convert the short
	// value to a double every time.
	m_dScale = _scale;
}

short CGeneticDataWriter::ConvertToShort(const double x) {
	short s1, s2;

	s1 = (unsigned short)floor(x * m_dScale);
	s2 = s1 + 1;
	if (fabs(x - ((double)s1)*m_dScale) > fabs(x - ((double)s2)*m_dScale))
		return s2;
	return s1;
}
////////////////////////////////////////////////////////////////////////////////
//							CDosageDataWriter
////////////////////////////////////////////////////////////////////////////////
// Class for writing genetic data to binary dosage files when only
// dosages need to be saved.

// Constructor
// Since only dosages are written the data written is always a vector of
// short integers of length the number of samples.
CDosageDataWriter::CDosageDataWriter(const unsigned short _scale, const unsigned int _sampleSize) : CGeneticDataWriter(_scale, _sampleSize) {
	m_dataToWrite.resize(m_sampleSize);
}

// Write the data to the file. Since only dosages are being writteng, the
// values of _p0, _p1, and _p2 are ignored. The data is always written at
// the end of the file.
int CDosageDataWriter::WriteData(std::fstream &_outfile, const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
	std::vector<unsigned short>::iterator usIt;
	std::vector<double>::const_iterator dIt;

	// Check if stream is capable of being written to.
	if (!_outfile.good())
		return 1;
	// Were the correct number of values passed?
	if (_dosage.size() != m_sampleSize)
		return 1;
	// Got ot the end of the file.
	_outfile.seekp(0, std::ios_base::end);

	usIt = m_dataToWrite.begin();
	// Loop over dosages
	for (dIt = _dosage.begin(); dIt != _dosage.end(); ++dIt, ++usIt) {
		if (*dIt < 0. || *dIt > 2.) {
			// Dosage value is not valid
			return 1;
		} else if (*dIt == *dIt) {
			*usIt = ConvertToShort(*dIt);
		} else {
			// Missing dosage value
			*usIt = 0xffff;
		}
	}
	// Write data
	_outfile.write((char *)m_dataToWrite.data(), m_sampleSize * sizeof(unsigned short));

	// Was data written successfully?
	if (!_outfile.good())
		return 1;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//							CGeneticDataWriter1
////////////////////////////////////////////////////////////////////////////////
// Class for writing genetic data to binary dosage files. Only the values for
// Pr(g=1) and Pr(g=2) are saved. P(g=0) and the dosage calculated from the
// saved values.

// Constructor
// Since only Pr(g=1) and Pr(g=2) are saved. The vector of values written has a
// length twice the number of samples.
CGeneticDataWriter1::CGeneticDataWriter1(const unsigned short _scale, const unsigned int _sampleSize) : CGeneticDataWriter(_scale, _sampleSize) {
	m_dataToWrite.resize(2 * m_sampleSize);
}

// Write the values to the file. Although only the values for _p1 and _p2 are
// written to the file, the values of _dosage and _p0 are used to check the
// validity fo the data. The values are written to the end of the file.
int CGeneticDataWriter1::WriteData(std::fstream &_outfile, const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
	std::vector<unsigned short>::iterator usIt1, usIt2;
	std::vector<double>::const_iterator dItd, dItp0, dItp1, dItp2;

	// Check if the file can be written to.
	if (!_outfile.good())
		return 1;
	// Are vectors for the appropriate sizes?
	if (_dosage.size() != m_sampleSize || _p0.size() != m_sampleSize || _p1.size() != m_sampleSize || _p2.size() != m_sampleSize)
		return 1;
	_outfile.seekp(0, std::ios_base::end);

	// Loop over the data.
	usIt1 = m_dataToWrite.begin();
	usIt2 = usIt1 + m_sampleSize;
	dItd = _dosage.begin();
	dItp0 = _p0.begin();
	dItp2 = _p2.begin();
	for (dItp1 = _p1.begin(); dItp1 != _p1.end(); ++dItd, ++dItp0, ++dItp1, ++dItp2, ++usIt1, ++usIt2) {
		// Is all the data not missing?
		if (*dItd == *dItd && *dItp0 == *dItp1 && *dItp0 == *dItp1 && *dItp2 == *dItp2) {
			// Are all the values in the required range?
			if (*dItd < 0. || *dItd > 2. || *dItp0 < 0. || *dItp0 > 1. || *dItp1 < 0. || *dItp1 > 1. || *dItp2 < 0. || *dItp2 > 1.)
				return 1;
			// Do the genetic probabilies add to 1 and does p1 + 2p2 = dosage?
			if (fabs(1. - *dItp0 - *dItp1 - *dItp2) > GeneticTolerance || fabs(*dItd - *dItp1 - *dItp2 - *dItp2) > GeneticTolerance)
				return 1;
			*usIt1 = ConvertToShort(*dItp1);
			*usIt2 = ConvertToShort(*dItp2);
		}
		// Is all the data missing?
		else if (!(*dItd == *dItd) && !(*dItp0 == *dItp0) && !(*dItp1 == *dItp1) && !(*dItp2 == *dItp2)) {
			// All missing
			*usIt1 = 0xffff;
			*usIt2 = 0xffff;
		}
		// Error - some data is missing and some is not
		else {
			return 1;
		}
	}
	// Write the data
	_outfile.write((char *)m_dataToWrite.data(), 2 * m_sampleSize * sizeof(unsigned short));
	// Was the data written successfully?
	if (!_outfile.good())
		return 1;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//							CGeneticDataWriter3
////////////////////////////////////////////////////////////////////////////////
// Class to write genetic data to binary dosage file where are the minimum
// number of values are stored such that all the values can be restored within
// the margin of error, 0.0001.

// Constructor
// Since it possible that all four values will be needed to keep required
// precision, the maximum size of the vector to be written out needs to be
// a vector of 4 times the number of samples.
CGeneticDataWriter3::CGeneticDataWriter3(const unsigned short _scale, const unsigned int _sampleSize) : CGeneticDataWriter(_scale, _sampleSize) {
	m_dataToWrite.resize(4 * m_sampleSize);
}

int CGeneticDataWriter3::WriteData(std::fstream &_outfile, const std::vector<double> &_dosage, const std::vector<double> &_p0, const std::vector<double> &_p1, const std::vector<double> &_p2) {
	std::vector<unsigned short>::iterator usIt1, usIt2;
	std::vector<double>::const_iterator dItd, dItp0, dItp1, dItp2;
	int numAdded;
	int outputLength;

	// Check if the file can be written to.
	if (!_outfile.good())
		return 1;
	// Are vectors for the appropriate sizes?
	if (_dosage.size() != m_sampleSize || _p0.size() != m_sampleSize || _p1.size() != m_sampleSize || _p2.size() != m_sampleSize)
		return 1;
	_outfile.seekp(0, std::ios_base::end);

	// Loop over the data.
	numAdded = 0;
	usIt1 = m_dataToWrite.begin();
	usIt2 = usIt1 + m_sampleSize;
	dItd = _dosage.begin();
	dItp0 = _p0.begin();
	dItp2 = _p2.begin();
	for (dItp1 = _p1.begin(); dItp1 != _p1.end(); ++dItd, ++dItp0, ++dItp1, ++dItp2, ++usIt1) {
		// Is all the data not missing?
		if (*dItd == *dItd && *dItp0 == *dItp0 && *dItp1 == *dItp1 && *dItp2 == *dItp2) {
			// Are all the values in the required range?
			if (*dItd < 0. || *dItd > 2. || *dItp0 < 0. || *dItp0 > 1. || *dItp1 < 0. || *dItp1 > 1. || *dItp2 < 0. || *dItp2 > 1.)
				return 1;
			// Do the genetic probabilies add to 1 and does p1 + 2p2 = dosage?
			if (fabs(1. - *dItp0 - *dItp1 - *dItp2) > GeneticTolerance || fabs(*dItd - *dItp1 - *dItp2 - *dItp2) > GeneticTolerance)
				return 1;
			*usIt1 = ConvertToShort(*dItd);
			// Can the values be derived within tolerances from a subset of values
			if (fabs(1. - *dItp0 - *dItp1 - *dItp2) > WritingTolerance || fabs(*dItd - *dItp1 - *dItp2 - *dItp2) > WritingTolerance) {
				// No, all four values must be written
				*usIt1 |= 0x8000;
				*usIt2 = (ConvertToShort(*dItp1) | 0x8000);
				++usIt2;
				*usIt2 = ConvertToShort(*dItp0);
				++usIt2;
				*usIt2 = ConvertToShort(*dItp2);
				++usIt2;
				numAdded += 3;
			}
			// Are Pr(g=0) and Pr(g=2) both non-zero?
			else if (*dItp0 != 0 && *dItp2 != 0) {
				// Yes, need to write dosage and p1
				*usIt1 |= 0x8000;
				*usIt2 = ConvertToShort(*dItp1);
				++usIt2;
				++numAdded;
			}
		}
		// Is all the data missing?
		else if (!(*dItd == *dItd) && !(*dItp0 == *dItp0) && !(*dItp1 == *dItp1) && !(*dItp2 == *dItp2)) {
			// All missing
			*usIt1 = 0xffff;
		}
		// Error - some data is missing and some is not
		else {
			return 1;
		}
	}
	// Calculate the length of the vector that needs to be writting (in bytes)
	// And write it to the file
	outputLength = (m_sampleSize + numAdded) * sizeof(short);
	_outfile.write((char *)&outputLength, sizeof(int));
	// Write the data
	_outfile.write((char *)m_dataToWrite.data(), outputLength);
	// Was the data written successfully?
	if (!_outfile.good())
		return 1;
	return 0;
}
