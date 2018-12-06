#ifndef VCFMINIREADER_H
#define VCFMINIREADER_H 1

#ifndef MINIREADER_H
#include "MiniReader.h"
#endif

class CVCFMiniReader : public CMiniReader {
protected:
public:
  CVCFMiniReader(const std::string &_filename, const int _numSamples, const int _numSNPs, const std::streampos _startData);
  virtual ~CVCFMiniReader() {}

  virtual int OpenFile();
};

#endif
