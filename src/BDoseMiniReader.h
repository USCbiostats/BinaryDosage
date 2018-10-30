#ifndef BDOSEMINIREADER_H
#define BDOSEMINIREADER_H 1

#ifndef MINIREADER_H
#include "MiniReader.h"
#endif

int GetBDoseFormat(const std::string &_filename, int &_version, int &_subversion);

class CBDoseMiniReader : public CMiniReader {
protected:
  int m_version, m_subversion;

  CBDoseMiniReader(const std::string &_filename);
public:
  virtual ~CBDoseMiniReader() {}

  virtual int OpenFile();
};

class CBDoseMiniReader1 : public CBDoseMiniReader {
protected:
public:
  CBDoseMiniReader1(const std::string &_filename, const int _numSamples, const int _numSNPs);
  virtual ~CBDoseMiniReader1() {}
};

class CBDoseMiniReader4 : public CBDoseMiniReader {
protected:
public:
  CBDoseMiniReader4(const std::string &_filename);
  virtual ~CBDoseMiniReader4() {}
};

#endif
