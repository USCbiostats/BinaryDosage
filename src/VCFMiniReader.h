#ifndef VCFMINIREADER_H
#define VCFMINIREADER_H 1

#ifndef MINIREADER_H
#include "MiniReader.h"
#endif

class CVCFMiniReader : CMiniReader {
protected:

public:
  CVCFMiniReader(const std::string &_filename, const unsigned int _numSamples, const unsigned int _numSNPs, const std::streampos _startData);
  virtual ~CVCFMiniReader() {}
};

#endif
