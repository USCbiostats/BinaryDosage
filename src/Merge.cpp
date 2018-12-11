#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include "GetBDoseInfo.h"
#include "BDoseMiniReader.h"
#include "VCFMiniReader.h"
#include "BDoseWriter.h"

CBDoseMiniReader *OpenBDoseMiniReader(const Rcpp::List &bdInfo) {
  std::string filename = bdInfo["filename"];
  int format = bdInfo["format"];
  int numSub = bdInfo["NumSamples"];
  int numSNPs = bdInfo["NumSNPs"];
  std::vector<int> index = Rcpp::as< std::vector<int> >(bdInfo["Indices"]);
  std::vector<int>::iterator intIt;
  std::vector<std::streampos>::iterator spIt;
  std::streampos snpIndex;
  std::vector<std::streampos> posIndex;

  CBDoseMiniReader *miniReader = NULL;

  if (format < 4)
    miniReader = new CBDoseMiniReader1(filename, numSub, numSNPs);
  else
    miniReader = new CBDoseMiniReader4(filename);

  snpIndex = 0;
  if (index.size() != 0) {
    posIndex.resize(index.size());
    for(intIt = index.begin(), spIt = posIndex.begin(); intIt != index.end(); ++intIt, ++spIt) {
      snpIndex += *intIt;
      *spIt = snpIndex;
    }
    miniReader->ChunkIt(posIndex);
  }

  return miniReader;
}

// [[Rcpp::export]]
int MergeBDC(const std::string &mergeFilename, const Rcpp::StringVector &filenames,
             const Rcpp::List &mergeInfo, const Rcpp::List &bdInfoList,
             const std::string &famFilename, const std::string &mapFilename,
             int format, int version) {
  int i;
  int numFiles;

  Rcpp::IntegerMatrix snpsToUse = mergeInfo["locations"];
  Rcpp::DataFrame subjects(mergeInfo["subjects"]);
  std::vector<std::string> fid = Rcpp::as<std::vector<std::string> >(subjects["fid"]);
  std::vector<std::string> sid = Rcpp::as<std::vector<std::string> >(subjects["iid"]);
  Rcpp::DataFrame snps(mergeInfo["snpsToMerge"]);
  std::vector<std::string> snpID = Rcpp::as<std::vector<std::string> >(snps["SNPID"]);
  std::vector<std::string> chromosome = Rcpp::as<std::vector<std::string> >(snps["Chromosome"]);
  std::vector<int> loc = Rcpp::as<std::vector<int> >(snps["Location"]);
  std::vector<std::string> ref = Rcpp::as<std::vector<std::string> >(snps["Reference"]);
  std::vector<std::string> alt = Rcpp::as<std::vector<std::string> >(snps["Alternate"]);
  std::vector<std::vector<double> > aaf, maf, avgCall, rSq;
  std::vector<std::vector<double> >::iterator vDblIt;
  std::vector<int> groupSizes;
  std::vector<double> dosage, p0, p1, p2;
  std::vector<double>::iterator itDosage, itP0, itP1, itP2;

  std::vector<CBDoseMiniReader *> bdmr;
  std::vector<CBDoseMiniReader *>::iterator bdmrIt;
  bool bdmrGood;
  CBDoseWriter *bdw;
  Rcpp::List x;
  int retVal = 0;

  // Create new file(s)
  if (format < 4)
    bdw = new CBDoseWriter1(mergeFilename, famFilename, mapFilename, format, version);
  else
    bdw = new CBDoseWriter4(mergeFilename, format, version);

  if (!bdw->good()) {
    Rcpp::Rcerr << "Unable to create new file" << std::endl;
    delete bdw;
  }

  // Open the files
  numFiles = bdInfoList.length();
  bdmr.resize(numFiles);
  for (i = 0; i < numFiles; ++i) {
    x = bdInfoList[i];
    bdmr[i] = OpenBDoseMiniReader(x);
    if (!(bdmr[i]->good())) {
      Rcpp::Rcerr << "Failed to open file\t" << i << std::endl;
      delete bdw;
      for (bdmrIt = bdmr.begin(); bdmrIt != bdmr.end(); ++bdmrIt) {
        if (*bdmrIt)
          delete *bdmrIt;
      }
    }
  }

  aaf.resize(chromosome.size());
  for (vDblIt = aaf.begin(); vDblIt != aaf.end(); ++vDblIt)
    vDblIt->resize(1);
  groupSizes.push_back(sid.size());

  bdmrGood = true;
  do {
    if (bdw->WriteGroupData(groupSizes)) {
      Rcpp::Rcerr << "Error writing group data" << std::endl;
      break;
    }
    if (bdw->WriteSubjectData(fid, sid)) {
      Rcpp::Rcerr << "Error writing subject data" << std::endl;
      break;
    }
    if (bdw->WriteSNPData(chromosome, snpID, loc, ref, alt, aaf, maf, avgCall, rSq)) {
      Rcpp::Rcerr << "Error writing SNP data" << std::endl;
      break;
    }

    dosage.resize(subjects.nrow());
    p0.resize(subjects.nrow());
    p1.resize(subjects.nrow());
    p2.resize(subjects.nrow());

    for (int i = 0; i < snpsToUse.nrow() && bdmrGood; ++i) {
      itDosage = dosage.begin();
      itP0 = p0.begin();
      itP1 = p1.begin();
      itP2 = p2.begin();
      for (int j = 0; j < snpsToUse.ncol(); ++j) {
        if (bdmr[j]->GetSNP(snpsToUse(i, j)) == false) {
          Rcpp::Rcerr << "Error reading SNP " << i << " from " << j << std::endl;
          bdmrGood = false;
          break;
        }
        for (int k = 0; k < bdmr[j]->NumSamples(); ++k, ++itDosage, ++itP0, ++itP1, ++itP2) {
          // Error-checking not really needed
          /*
          if (itDosage == dosage.end())
            Rcpp::Rcout << "Dosage pointer to big" << std::endl;
          if (itP0 == p0.end())
            Rcpp::Rcout << "P0 pointer to big" << std::endl;
          if (itP1 == p1.end())
            Rcpp::Rcout << "P1 pointer to big" << std::endl;
          if (itP2 == p2.end())
            Rcpp::Rcout << "P2 pointer to big" << std::endl;
          */
          *itDosage = bdmr[j]->Dosage()[k];
          *itP0 = bdmr[j]->P0()[k];
          *itP1 = bdmr[j]->P1()[k];
          *itP2 = bdmr[j]->P2()[k];
        }
      }
      if (bdmrGood) {
        if (bdw->WriteGeneticData(dosage, p0, p1, p2)) {
          Rcpp::Rcerr << "Error writing genetic data" << std::endl;
          bdmrGood = false;
        }
        else {
          aaf[i][0] = std::accumulate(dosage.begin(), dosage.end(), 0.) / (subjects.nrow() + subjects.nrow());
        }
      }
    }

    if (bdw->UpdateSNPInfo(aaf, maf, avgCall, rSq))
      bdmrGood = false;
  } while (0);
  if (bdmrGood) {
    if (bdw->Finalize())
      Rcpp::Rcerr << "Error finalizing file" << std::endl;
  } else {
    retVal = 1;
  }

  if (bdw)
    delete bdw;
  for (bdmrIt = bdmr.begin(); bdmrIt != bdmr.end(); ++bdmrIt) {
    if (*bdmrIt)
      delete *bdmrIt;
  }
  return retVal;
}

