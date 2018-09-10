#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
#include "VCFtoBinaryDosage.h"

CVCFtoBinaryDosage42::CVCFtoBinaryDosage42() {
  m_bd = NULL;
  m_vcfFile = NULL;
}

CVCFtoBinaryDosage42::~CVCFtoBinaryDosage42() {
  if (m_bd)
    delete [] m_bd;
  if (m_vcfFile)
    delete m_vcfFile;
}

void CVCFtoBinaryDosage42::OpenTempFiles() {
  std::string outFilename;

  outFilename = m_bdFilename + std::string(".tmp.chromosome");
  m_chromosomeFile.open(outFilename.c_str());
  outFilename = m_bdFilename + std::string(".tmp.bp");
  m_bpFile.open(outFilename.c_str(), std::ios_base::out | std::ios_base::binary);
  outFilename = m_bdFilename + std::string(".tmp.snpID");
  m_snpIDFile.open(outFilename.c_str());
  outFilename = m_bdFilename + std::string(".tmp.ref");
  m_refFile.open(outFilename.c_str());
  outFilename = m_bdFilename + std::string(".tmp.alt");
  m_altFile.open(outFilename.c_str());
  outFilename = m_bdFilename + std::string(".tmp.altFreq");
  m_altFreqFile.open(outFilename.c_str(), std::ios_base::out | std::ios_base::binary);
  outFilename = m_bdFilename + std::string(".tmp.bd");
  m_doseFile.open(outFilename.c_str(), std::ios_base::out | std::ios_base::binary);
}

void CVCFtoBinaryDosage42::CloseTempFiles() {
  m_chromosomeFile.close();
  m_bpFile.close();
  m_snpIDFile.close();
  m_refFile.close();
  m_altFile.close();
  m_altFreqFile.close();
  m_doseFile.close();
}

void CVCFtoBinaryDosage42::WriteDosage(const std::vector<double> &d, const std::vector<double> &p0,
                                         const std::vector<double> &p1, const std::vector<double> &p2) {
  int numSub, numAdded;
  const double *dp, *pp0, *pp1, *pp2;
  short sd, sp0, sp1, sp2;
  unsigned short *pbd, *pbdex;
  int i;

  numSub = d.size();
  dp = &d[0];
  pp0 = &p0[0];
  pp1 = &p1[0];
  pp2 = &p2[0];
  pbd = m_bd;
  pbdex = pbd + numSub;
  numAdded = 0;

  for (i = 0; i < numSub; ++i, ++dp, ++pp0, ++pp1, ++pp2, ++pbd) {
    sd = 10000 * *dp;
//    if (i < 5 || i == numSub - 1)
//      Rcpp::Rcout << fabs(*dp - sd/10000.) << '\t' << fabs(*dp - (sd + 1)/10000.) << '\t';
    if (fabs(*dp - sd/10000.) > fabs(*dp - (sd + 1)/10000.))
      ++sd;
    sp0 = 10000 * *pp0;
    if (fabs(*pp0 - sp0/10000.) > fabs(*pp0 - (sp0 + 1)/10000.))
      ++sp0;
    sp1 = 10000 * *pp1;
    if (fabs(*pp1 - sp1/10000.) > fabs(*pp1 - (sp1 + 1)/10000.))
      ++sp1;
    sp2 = 10000 * *pp2;
    if (fabs(*pp2 - sp2/10000.) > fabs(*pp2 - (sp2 + 1)/10000.))
      ++sp2;
    *pbd = sd;
    if (sp1 + sp2 + sp2 != sd || sp0 + sp1 + sp2 != 10000) {
      *pbd |= 0x8000;
      *pbdex = sp1;
      *pbdex |= 0x8000;
      ++pbdex;
      *pbdex = sp0;
      ++pbdex;
      *pbdex = sp2;
      ++pbdex;
      numAdded += 3;
    } else if (sp0 != 0 && sp2 != 0) {
      *pbd |= 0x8000;
      *pbdex = sp1;
      ++pbdex;
      numAdded += 1;
    }
//    if (i < 5 || i == numSub - 1)
//      Rcpp::Rcout << *dp << '\t' << sd/10000. << '\t' << sp0/10000. << '\t' << sp1/10000. << '\t' << sp2/10000. << std::endl;
  }
//  Rcpp::Rcout << std::endl;
  numAdded += numSub;
  numAdded *= sizeof(unsigned short);
  m_doseFile.write((char *)&numAdded, sizeof(int));
  m_doseFile.write((char *)m_bd, numAdded);
}

void CVCFtoBinaryDosage42::WriteHeader() {
  const char header[8] = { 'b', 'o', 's', 'e', 0x0, 0x4, 0x0, 0x2};
  const int zero = 0;
  const int numGroups = 1;

  m_binaryDosageFile.write(header, 8);
  m_binaryDosageFile.write((char *)&m_numSub, sizeof(int));
  m_binaryDosageFile.write((char *)&m_numSNPs, sizeof(int));
  m_binaryDosageFile.write((char *)&numGroups, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));

}

void CVCFtoBinaryDosage42::WriteGroups() {
  std::streampos startSubjects;
  int startSub;

  m_binaryDosageFile.write((char *)&m_numSub, sizeof(int));
  startSubjects = m_binaryDosageFile.tellp();
  startSub = startSubjects;
  m_binaryDosageFile.seekp(28);
  m_binaryDosageFile.write((char *)&startSub, sizeof(int));
  m_binaryDosageFile.seekp(startSubjects);
}

void CVCFtoBinaryDosage42::WriteSubjects(const std::vector<std::string> &subID) {
  const int zero = 0;
  const char zeroCh = 0;
  std::streampos startSub, startSNP;
  int startSNP2;
  std::ostringstream oss;
  int subSize;

  startSub = m_binaryDosageFile.tellp();
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  m_binaryDosageFile.write((char *)&zero, sizeof(int));
  for (std::vector<std::string>::const_iterator sub = subID.begin(); sub != subID.end(); ++sub)
    oss << *sub << '\t';
  m_binaryDosageFile.write((char *)oss.str().c_str(), oss.str().size() - 1);
  m_binaryDosageFile.write(&zeroCh, 1);
  subSize = oss.str().size();
  startSNP = m_binaryDosageFile.tellp();
  m_binaryDosageFile.seekp(startSub);
  m_binaryDosageFile.write((char *)&subSize, sizeof(int));
  m_binaryDosageFile.seekp(startSNP);
  startSNP2 = startSNP;
  m_binaryDosageFile.seekp(32);
  m_binaryDosageFile.write((char *)&startSNP2, sizeof(int));
  m_binaryDosageFile.seekp(startSNP);
}

void CVCFtoBinaryDosage42::WriteSNPInfo(const int snpOptions, const std::string &singleChromosome) {
  int nameSize = 0;
  int chromosomeSize = 0;
  int refSize = 0;
  int altSize = 0;
  char zero = 0x0;
  std::ifstream infile;
  std::string filename;
  std::string fileContents;
  int *bp = NULL;
  double *dval = NULL;
  std::streampos startSNP, startDosage, fileSize;
  int startDosage2;

  startSNP = m_binaryDosageFile.tellp();
  m_binaryDosageFile.write((char *)&nameSize, sizeof(int));
  m_binaryDosageFile.write((char *)&chromosomeSize, sizeof(int));
  m_binaryDosageFile.write((char *)&refSize, sizeof(int));
  m_binaryDosageFile.write((char *)&altSize, sizeof(int));
  if (snpOptions & 0x0002) {
    filename = m_bdFilename + std::string(".tmp.snpID");
    infile.open(filename.c_str());
    getline(infile, fileContents);
    m_binaryDosageFile.write(fileContents.c_str(), fileContents.size() - 1);
    m_binaryDosageFile.write(&zero, 1);
    nameSize = fileContents.size();
    infile.close();
    infile.clear();
  }
  if (snpOptions & 0x0004) {
    if (snpOptions & 0x0008) {
      m_binaryDosageFile.write(singleChromosome.c_str(), singleChromosome.size());
      m_binaryDosageFile.write(&zero, 1);
      chromosomeSize = singleChromosome.size() + 1;
    } else {
      filename = m_bdFilename + std::string(".tmp.chromosome");
      infile.open(filename.c_str());
      getline(infile, fileContents);
      m_binaryDosageFile.write(fileContents.c_str(), fileContents.size() - 1);
      m_binaryDosageFile.write(&zero, 1);
      chromosomeSize = fileContents.size();
      infile.close();
      infile.clear();
    }
  }
  if (snpOptions & 0x0010) {
    bp = new int[m_numSNPs];
    filename = m_bdFilename + std::string(".tmp.bp");
    infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    infile.read((char *)bp, m_numSNPs * sizeof(int));
    m_binaryDosageFile.write((char *)bp, m_numSNPs * sizeof(int));
    infile.close();
    infile.clear();
  }
  if (snpOptions & 0x0020) {
    filename = m_bdFilename + std::string(".tmp.ref");
    infile.open(filename.c_str());
    if (!infile.good()) {
      Rcpp::Rcerr << "Unable to open file " << filename << std::endl;
      return;
    }
    getline(infile, fileContents);
    m_binaryDosageFile.write(fileContents.c_str(), fileContents.size() - 1);
    m_binaryDosageFile.write(&zero, 1);
    refSize = fileContents.size();
    infile.close();
    infile.clear();
  }
  if (snpOptions & 0x0040) {
    filename = m_bdFilename + std::string(".tmp.alt");
    infile.open(filename.c_str());
    getline(infile, fileContents);
    m_binaryDosageFile.write(fileContents.c_str(), fileContents.size() - 1);
    m_binaryDosageFile.write(&zero, 1);
    altSize = fileContents.size();
    infile.close();
  }
  if (snpOptions & 0x0080) {
    if (dval == NULL)
      dval = new double[m_numSNPs];
    filename = m_bdFilename + std::string(".tmp.altFreq");
    infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    infile.read((char *)dval, m_numSNPs * sizeof(double));
    m_binaryDosageFile.write((char *)dval, m_numSNPs * sizeof(double));
    infile.close();
  }
  if (snpOptions & 0x0100) {
    if (dval == NULL)
      dval = new double[m_numSNPs];
    filename = m_bdFilename + std::string(".tmp.maf");
    infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    infile.read((char *)dval, m_numSNPs * sizeof(double));
    m_binaryDosageFile.write((char *)dval, m_numSNPs * sizeof(double));
    infile.close();
  }
  if (snpOptions & 0x0200) {
    if (dval == NULL)
      dval = new double[m_numSNPs];
    filename = m_bdFilename + std::string(".tmp.avg");
    infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    infile.read((char *)dval, m_numSNPs * sizeof(double));
    m_binaryDosageFile.write((char *)dval, m_numSNPs * sizeof(double));
    infile.close();
  }
  if (snpOptions & 0x0400) {
    if (dval == NULL)
      dval = new double[m_numSNPs];
    filename = m_bdFilename + std::string(".tmp.rsq");
    infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    infile.read((char *)dval, m_numSNPs * sizeof(double));
    m_binaryDosageFile.write((char *)dval, m_numSNPs * sizeof(double));
    infile.close();
  }
  if (bp)
    delete [] bp;
  if (dval)
    delete [] dval;

  startDosage = m_binaryDosageFile.tellp();
  m_binaryDosageFile.seekp(startSNP);
  m_binaryDosageFile.write((char *)&nameSize, sizeof(int));
  m_binaryDosageFile.write((char *)&chromosomeSize, sizeof(int));
  m_binaryDosageFile.write((char *)&refSize, sizeof(int));
  m_binaryDosageFile.write((char *)&altSize, sizeof(int));
  m_binaryDosageFile.seekp(24);
  m_binaryDosageFile.write((char *)&snpOptions, sizeof(int));
  m_binaryDosageFile.seekp(36);
  startDosage2 = startDosage;
  m_binaryDosageFile.write((char *)&startDosage2, sizeof(int));
  m_binaryDosageFile.seekp(startDosage);
}

void CVCFtoBinaryDosage42::WriteGeneticValues() {
  int snpSize;
  int i;
  std::ifstream infile;
  std::string filename;

  filename = m_bdFilename + std::string(".tmp.bd");
  infile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  for (i = 0; i < m_numSNPs; ++i) {
    infile.read((char *)&snpSize, sizeof(int));
    m_binaryDosageFile.write((char *)&snpSize, sizeof(int));
    infile.read((char *)m_bd, snpSize);
    m_binaryDosageFile.write((char *)m_bd, snpSize);
  }
  infile.close();
}

int CVCFtoBinaryDosage42::OpenVCF(const std::string &vcfFilename) {
  if (m_vcfFile)
    delete m_vcfFile;
  m_vcfFile = NULL;
  if (m_bd)
    delete [] m_bd;
  m_bd = NULL;

  m_vcfFile = new CReadVCF_Generic(vcfFilename);

  return 0;
}

int CVCFtoBinaryDosage42::Convert(const std::string &vcfFilename, const std::string &bdFilename) {
  bool oneChromosome;
  std::string singleChromosome;
  double altFreq;

  m_bdFilename = bdFilename;
  OpenVCF(vcfFilename);
  if (m_vcfFile->ReadSubjects()) {
    Rcpp::Rcerr << "VCF read error" << std::endl;
    return 1;
  }

  m_numSub = m_vcfFile->NumSubjects();
  m_bd = new unsigned short[4 * m_numSub];

  if (m_vcfFile->GetFirstSNP()) {
    Rcpp::Rcout << "Error reading first SNP" << std::endl;
    return 1;
  }
  oneChromosome = true;
  singleChromosome = "";
  OpenTempFiles();
  do {
    if (singleChromosome == "")
      singleChromosome = m_vcfFile->Chromosome()[0];
    else if (singleChromosome != m_vcfFile->Chromosome()[0])
      oneChromosome = false;
    m_chromosomeFile << m_vcfFile->Chromosome()[0] << '\t';
    m_bpFile.write((char *)&m_vcfFile->Location()[0], sizeof(int));
    m_snpIDFile << m_vcfFile->SNPID()[0] << '\t';
    m_refFile << m_vcfFile->ReferenceAllele()[0] << '\t';
    m_altFile << m_vcfFile->AlternateAllele()[0] << '\t';
    altFreq = m_vcfFile->AlternateAlleleFrequency();
    m_altFreqFile.write((char *)&altFreq, sizeof(double));
    WriteDosage(m_vcfFile->Dosage(), m_vcfFile->P0(), m_vcfFile->P1(), m_vcfFile->P2());
  } while (!m_vcfFile->GetNextSNP());
  m_numSNPs = m_vcfFile->NumSNPs();
  CloseTempFiles();

  m_binaryDosageFile.open(m_bdFilename.c_str(), std::ios_base::out | std::ios_base::binary);
  WriteHeader();
  WriteGroups();
  WriteSubjects(m_vcfFile->SubjectID());
  if (oneChromosome)
    WriteSNPInfo(0x00fc, singleChromosome);
  else
    WriteSNPInfo(0x00f4, singleChromosome);
  WriteGeneticValues();
  m_binaryDosageFile.close();
//  OpenTempFiles();
//  CloseTempFiles();

  return 0;
}

int CVCF53toBinaryDosage42::OpenVCF(const std::string &vcfFilename) {
  if (m_vcfFile)
    delete m_vcfFile;
  m_vcfFile = NULL;
  if (m_bd)
    delete [] m_bd;
  m_bd = NULL;

  m_vcfFile = new CReadVCF_HRC(vcfFilename);

  return 0;
}

void CVCF53toBinaryDosage42::WriteDosage(const std::vector<double> &d, const std::vector<double> &p0,
                                         const std::vector<double> &p1, const std::vector<double> &p2) {
  int numSub, numAdded;
  const double *dp, *pp0, *pp1, *pp2;
  unsigned short *pbd, *pbdex;
  int i;

  numSub = d.size();
  dp = &d[0];
  pp0 = &p0[0];
  pp1 = &p1[0];
  pp2 = &p2[0];
  pbd = m_bd;
  pbdex = pbd + numSub;
  numAdded = 0;

  for (i = 0; i < numSub; ++i, ++dp, ++pp0, ++pp1, ++pp2, ++pbd) {
    *pbd = 10000 * (*dp + 1e-6);
    if (fabs(*dp - (*pp1 + *pp2 + *pp2)) > 1e-6 || fabs(1. - (*pp0 + *pp1 + *pp2)) > 1e-6) {
      *pbd |= 0x8000;
      *pbdex = 10000 * (*pp1 + 1e-6);
      *pbdex |= 0x8000;
      ++pbdex;
      *pbdex = 10000 * (*pp0 + 1e-6);
      ++pbdex;
      *pbdex = 10000 * (*pp2 + 1e-6);
      ++pbdex;
      numAdded += 3;
    } else if (fabs(*pp0) > 1e-6 && fabs(*pp2) > 1e-6) {
      *pbd |= 0x8000;
      *pbdex = 10000 * (*pp1 + 1e-6);
      ++pbdex;
      numAdded += 1;
    }
  }
  numAdded += numSub;
  numAdded *= sizeof(unsigned short);
  m_doseFile.write((char *)&numAdded, sizeof(int));
  m_doseFile.write((char *)m_bd, numAdded);
}
