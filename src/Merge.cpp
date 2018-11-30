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
  CBDoseMiniReader *bdmr = NULL;

  if (format < 4)
    bdmr = new CBDoseMiniReader1(filename, numSub, numSNPs);
  else
    bdmr = new CBDoseMiniReader4(filename);

  return bdmr;
}

// [[Rcpp::export]]
int MergeBDC(const std::string &mergeFilename, const Rcpp::StringVector &filenames,
             const Rcpp::List &mergeInfo, const Rcpp::List &bdInfoList,
             const std::string &famFilename, const std::string &mapFilename,
             int format, int version, int batchSize) {
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
  std::vector<double>::iterator dblIt1, dblIt2, dblIt3, dblIt4;
  std::vector<int> groupSizes;
  std::vector<std::vector<double> > dosage, p0, p1, p2;
  std::vector<std::vector<double> >::iterator itDosage, itP0, itP1, itP2;

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
      Rcpp::Rcout << "Failed to open file\t" << i << std::endl;
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

    dosage.resize(batchSize);
    p0.resize(batchSize);
    p1.resize(batchSize);
    p2.resize(batchSize);
    itDosage = dosage.begin();
    itP0 = p0.begin();
    itP1 = p1.begin();
    itP2 = p2.begin();
    for(;itDosage != dosage.end(); ++itDosage, ++itP0, ++itP1, ++itP2) {
      itDosage->resize(subjects.nrow());
      itP0->resize(subjects.nrow());
      itP1->resize(subjects.nrow());
      itP2->resize(subjects.nrow());
    }
    itDosage = dosage.begin();
    itP0 = p0.begin();
    itP1 = p1.begin();
    itP2 = p2.begin();

    for (int i = 0; i < snpsToUse.nrow() && bdmrGood; ++i) {
      Rcpp::Rcout << i << '\t';
      dblIt1 = itDosage->begin();
      dblIt2 = itP0->begin();
      dblIt3 = itP1->begin();
      dblIt4 = itP2->begin();

      for (int j = 0; j < snpsToUse.ncol(); ++j) {
        if (bdmr[j]->GetSNP(snpsToUse(i, j)) == false) {
          Rcpp::Rcout << "Error reading SNP " << i << " from " << j << std::endl;
          bdmrGood = false;
          break;
        }
        for (int k = 0; k < bdmr[j]->NumSamples(); ++k, ++dblIt1, ++dblIt2, ++dblIt3, ++dblIt4) {
          *dblIt1 = bdmr[j]->Dosage()[k];
          *dblIt2 = bdmr[j]->P0()[k];
          *dblIt3 = bdmr[j]->P1()[k];
          *dblIt4 = bdmr[j]->P2()[k];
        }
        if (bdmrGood) {
          Rcpp::Rcout << dosage[i][0] << '\t' << p0[i][0] << '\t' << p1[i][0] << '\t' << p2[i][0] << std::endl;
          if (bdw->WriteGeneticData(dosage[i], p0[i], p1[i], p2[i])) {
            Rcpp::Rcout << "Error writing genetic data" << std::endl;
            bdmrGood = false;
          }
          else {
            aaf[i][0] = std::accumulate(dosage[i].begin(), dosage[i].end(), 0.) / (subjects.nrow() + subjects.nrow());
          }
        }
      }
    }
    Rcpp::Rcout << std::endl;

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

