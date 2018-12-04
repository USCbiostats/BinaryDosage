#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include "GeneticDataReader.h"

////////////////////////////////////////////////////////////////////////////////
//							CBDoseDataReader
////////////////////////////////////////////////////////////////////////////////
// Base class for reading genetic data from binary dosage files
// When it is created the scale parameter and number of observations
// must be specified. These are constants and cannot be changed.

// Constructor
CBDoseDataReader::CBDoseDataReader(const unsigned short _scale, const int _sampleSize) : CGeneticDataReader(_sampleSize), m_scale(_scale) {
  // The scale is changed to a double. This helps the speed of
  // execution by not having to constantly convert the short
  // value to a double every time.
  m_dScale = _scale;
}

////////////////////////////////////////////////////////////////////////////////
//							CBDoseDosageReader
////////////////////////////////////////////////////////////////////////////////
// Class for reading genetic data from binary dosage files when only
// dosages are saved.

// Constructor
// Since only dosages are read, the data read is always a vector of
// short integers of length the number of samples.
CBDoseDosageReader::CBDoseDosageReader(const unsigned short _scale, const int _sampleSize) : CBDoseDataReader(_scale, _sampleSize) {
  m_dataToRead.resize(m_sampleSize);
}

// Read the data to the file. Since only dosages are being read, the
// values of _p0, _p1, and _p2 are ignored. The data is always read from
// the current file location.
int CBDoseDosageReader::ReadData(std::ifstream &_infile, std::vector<double> &_dosage, std::vector<double> &_p0, std::vector<double> &_p1, std::vector<double> &_p2) {
  std::vector<unsigned short>::const_iterator usIt;
  std::vector<double>::iterator dIt;

  // Check if stream is capable of being read from?
  if (!_infile.good())
    return 1;
  // Were the correct number of values passed?
  if ((int)_dosage.size() != m_sampleSize)
    return 1;

  // Read in the data
  _infile.read((char *)m_dataToRead.data(), m_sampleSize * sizeof(unsigned short));
  // Did data read in successfully?
  if (!_infile.good())
    return 1;

  usIt = m_dataToRead.begin();
  // Loop over dosages
  for (dIt = _dosage.begin(); dIt != _dosage.end(); ++dIt, ++usIt) {
    if (*usIt == 0xffff)
      // Missing
      *dIt = std::numeric_limits<double>::quiet_NaN();
    else
      *dIt = *usIt / m_dScale;
  }

  return 0;
}

// Skip over the next SNP
int CBDoseDosageReader::SkipSNP(std::ifstream &_infile) {
  // Check if stream is capable of being read from
  if (!_infile.good())
    return 1;
  // Skip the data for the next SNP
  _infile.seekg(m_sampleSize * sizeof(unsigned short), std::ios_base::cur);
  // Check if stream is still good
  if (!_infile.good())
    return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//							CBDose1DataReader
////////////////////////////////////////////////////////////////////////////////
// Class for reading genetic data from a binary dosage files. Only the values for
// Pr(g=1) and Pr(g=2) are saved. P(g=0) and the dosage calculated from the
// saved values.

// Constructor
// Since only Pr(g=1) and Pr(g=2) are saved. The vector of values read has a
// length twice the number of samples.
CBDose1DataReader::CBDose1DataReader(const unsigned short _scale, const int _sampleSize) : CBDoseDataReader(_scale, _sampleSize) {
  m_dataToRead.resize(2 * m_sampleSize);
}

// Write the values to the file. Although only the values for _p1 and _p2 are
// written to the file, the values of _dosage and _p0 are used to check the
// validity fo the data. The values are written to the end of the file.
int CBDose1DataReader::ReadData(std::ifstream &_infile, std::vector<double> &_dosage, std::vector<double> &_p0, std::vector<double> &_p1, std::vector<double> &_p2) {
  std::vector<unsigned short>::const_iterator usIt1, usIt2;
  std::vector<double>::iterator dItd, dItp0, dItp1, dItp2;

  // Check if the file can be read from.
  if (!_infile.good())
    return 1;
  // Are vectors for the appropriate sizes?
  if ((int)_dosage.size() != m_sampleSize || (int)_p0.size() != m_sampleSize || (int)_p1.size() != m_sampleSize || (int)_p2.size() != m_sampleSize)
    return 1;

  // Write the data
  _infile.read((char *)m_dataToRead.data(), 2 * m_sampleSize * sizeof(unsigned short));
  // Was the data written successfully?
  if (!_infile.good())
    return 1;

  // Loop over the data.
  usIt1 = m_dataToRead.begin();
  usIt2 = usIt1 + m_sampleSize;
  dItd = _dosage.begin();
  dItp0 = _p0.begin();
  dItp2 = _p2.begin();
  for (dItp1 = _p1.begin(); dItp1 != _p1.end(); ++dItd, ++dItp0, ++dItp1, ++dItp2, ++usIt1, ++usIt2) {
    if (*usIt1 == 0xffff) {
      // Missing
      *dItd = std::numeric_limits<double>::quiet_NaN();
      *dItp0 = std::numeric_limits<double>::quiet_NaN();
      *dItp1 = std::numeric_limits<double>::quiet_NaN();
      *dItp2 = std::numeric_limits<double>::quiet_NaN();
    }
    else {
      *dItp1 = *usIt1 / m_dScale;
      *dItp2 = *usIt2 / m_dScale;
      *dItp0 = 1. - *dItp1 - *dItp2;
      if (*dItp0 < 0.)
        *dItp0 = 0.;
      *dItd = *dItp1 + *dItp2 + *dItp2;
      if (*dItd > 2.)
        *dItd = 2.;
    }
  }
  return 0;
}

// Skip over the next SNP
int CBDose1DataReader::SkipSNP(std::ifstream &_infile) {
  // Check if stream is capable of being read from
  if (!_infile.good())
    return 1;
  // Skip the data for the next SNP
  _infile.seekg(2 * m_sampleSize * sizeof(unsigned short), std::ios_base::cur);
  // Check if stream is still good
  if (!_infile.good())
    return 1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//							CBDose3DataReader
////////////////////////////////////////////////////////////////////////////////
// Class to read genetic data from a binary dosage file where the minimum
// number of values are stored such that all the values can be restored within
// the margin of error, 0.0001.

// Constructor
// Since it possible that all four values will be needed to keep required
// precision, the maximum size of the vector to be read needs to be
// a vector of 4 times the number of samples.
CBDose3DataReader::CBDose3DataReader(const unsigned short _scale, const int _sampleSize) : CBDoseDataReader(_scale, _sampleSize) {
  m_dataToRead.resize(4 * m_sampleSize);
}

int CBDose3DataReader::ReadData(std::ifstream &_infile, std::vector<double> &_dosage, std::vector<double> &_p0, std::vector<double> &_p1, std::vector<double> &_p2) {
  std::vector<unsigned short>::const_iterator usIt1, usIt2;
  std::vector<double>::iterator dItd, dItp0, dItp1, dItp2;
  int inputLength;

  // Check if the file can be read from.
  if (!_infile.good())
    return 1;
  // Are vectors for the appropriate sizes?
  if ((int)_dosage.size() != m_sampleSize || (int)_p0.size() != m_sampleSize || (int)_p1.size() != m_sampleSize || (int)_p2.size() != m_sampleSize)
    return 1;

  _infile.read((char *)&inputLength, sizeof(int));
  // Write the data
  _infile.read((char *)m_dataToRead.data(), inputLength);
  // Was the data written successfully?
  if (!_infile.good())
    return 1;

  // Loop over the data.
  usIt1 = m_dataToRead.begin();
  usIt2 = usIt1 + m_sampleSize;
  dItd = _dosage.begin();
  dItp0 = _p0.begin();
  dItp2 = _p2.begin();
  for (dItp1 = _p1.begin(); dItp1 != _p1.end(); ++dItd, ++dItp0, ++dItp1, ++dItp2, ++usIt1) {
    if (*usIt1 == 0xffff) {
      // Missing
      *dItd = std::numeric_limits<double>::quiet_NaN();
      *dItp0 = std::numeric_limits<double>::quiet_NaN();
      *dItp1 = std::numeric_limits<double>::quiet_NaN();
      *dItp2 = std::numeric_limits<double>::quiet_NaN();
    }
    else {
      if (*usIt1 & 0x8000) {
        // At least two values written
        *dItd = (*usIt1 & 0x7fff) / m_dScale;
        if (*usIt2 & 0x8000) {
          // Four values written
          *dItp1 = (*usIt2 & 0x7fff) / m_dScale;
          ++usIt2;
          *dItp0 = *usIt2 / m_dScale;
          ++usIt2;
          *dItp2 = *usIt2 / m_dScale;
          ++usIt2;
        }
        else {
          // Two values written
          *dItp1 = *usIt2 / m_dScale;
          ++usIt2;
          *dItp2 = (*dItd - *dItp1) / 2.;
          if (*dItp2 < 0.)
            *dItp2 = 0.;
          *dItp0 = 1. - *dItp1 - *dItp2;
          if (*dItp0 < 0.)
            *dItp0 = 0.;
        }
      }
      else {
        // Only dosage written
        *dItd = *usIt1 / m_dScale;
        if (*dItd < 1.) {
          *dItp1 = *dItd;
          *dItp0 = 1. - *dItp1;
          *dItp2 = 0.;
        }
        else {
          *dItp2 = *dItd - 1.;
          *dItp1 = 1. - *dItp2;
          *dItp0 = 0.;
        }
      }
    }
  }

  return 0;
}

// Skip over the next SNP
int CBDose3DataReader::SkipSNP(std::ifstream &_infile) {
  int snpSize;
  // Check if stream is capable of being read from
  if (!_infile.good())
    return 1;
  _infile.read((char *)&snpSize, sizeof(int));
  // Skip the data for the next SNP
  _infile.seekg(snpSize, std::ios_base::cur);
  // Check if stream is still good
  if (!_infile.good())
    return 1;
  return 0;
}

int CVCFDataReader::ReadData(std::ifstream &_infile, std::vector<double> &_dosage, std::vector<double> &_p0, std::vector<double> &_p1, std::vector<double> &_p2) {
  std::istringstream iss, iss2;
  std::string format;
  std::string formatVal;
  std::string dataStr;
  std::string values;
  double p0, p1, p2;
  int numFormat;
  int i;
  int dsLoc, gpLoc;
  std::vector<double>::iterator dIt, p0It, p1It, p2It;

  _infile >> format >> format >> format >> format >> format >> format >> format >> format >> format;
  std::replace(format.begin(), format.end(), ':', ' ');

  dsLoc = -1;
  gpLoc = -1;
  numFormat = 0;
  iss.str(format);
  iss >> formatVal;
  do {
    if (formatVal == "GP")
      gpLoc = numFormat;
    if (formatVal == "DS")
      dsLoc = numFormat;
    iss >> formatVal;
    ++numFormat;
  } while (!iss.fail());

  if (dsLoc == -1) {
    std::fill(_dosage.begin(), _dosage.end(), std::numeric_limits<double>::quiet_NaN());
  }
  if (gpLoc == -1) {
    std::fill(_p0.begin(), _p0.end(), std::numeric_limits<double>::quiet_NaN());
    std::fill(_p1.begin(), _p1.end(), std::numeric_limits<double>::quiet_NaN());
    std::fill(_p2.begin(), _p2.end(), std::numeric_limits<double>::quiet_NaN());
  }

  p0It = _p0.begin();
  p1It = _p1.begin();
  p2It = _p2.begin();
  for (dIt = _dosage.begin(); dIt != _dosage.end(); ++dIt, ++p0It, ++p1It, ++p2It) {
    _infile >> dataStr;
    std::replace(dataStr.begin(), dataStr.end(), ':', ' ');
    iss.clear();
    iss.str(dataStr);
    for (i = 0; i < numFormat; ++i) {
      iss >> values;
      if (i == dsLoc) {
        *dIt = std::stod(values);
        if (*dIt > 2.01 || *dIt < -0.01)
          return 1;
        if (*dIt > 2.)
          *dIt = 2.;
        else if (*dIt < 0.)
          *dIt = 0;
      }
      if (i == gpLoc) {
        std::replace(values.begin(), values.end(), ',', ' ');
        iss2.clear();
        iss2.str(values);
        iss2 >> p0 >> p1 >> p2;
        *p0It = p0;
        *p1It = p1;
        *p2It = p2;
        if (*p0It > 1.01 || *p0It < -0.01 || *p1It > 1.01 || *p1It < -0.01 || *p2It > 1.01 || *p2It < -0.01)
          return 1;
        if (*p0It > 1.)
          *p0It = 1.;
        else if (*p0It < 0.)
          *p0It = 0.;
        if (*p1It > 1.)
          *p1It = 1.;
        else if (*p1It < 0.)
          *p1It = 0.;
        if (*p2It > 1.)
          *p2It = 1.;
        else if (*p2It < 0.)
          *p2It = 0.;
        if(fabs(1. - *p0It - *p1It - *p2It) > 0.01)
          return 1;
      }
    }
    if (gpLoc != -1 && dsLoc != -1) {
      if (fabs(*dIt - *p1It - *p2It - *p2It) > 0.01)
        return 1;
    }
  }

  return 0;
}

int CVCFDataReader::SkipSNP(std::ifstream &_infile) {
  std::string junk;

  std::getline(_infile, junk);
  if (!_infile.good())
    return 1;
  return 0;
}


