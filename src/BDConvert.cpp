#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <limits>
#include <algorithm>
#include "GetBDoseInfo.h"
#include "BDoseMiniReader.h"
#include "VCFMiniReader.h"
#include "BDoseWriter.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List BDConvertC(const Rcpp::List &bdInfo, const std::string &newFile, const std::string &famFile,
                      const std::string &mapFile, int newFormat, int newVersion) {
  Rcpp::List retVal;
  std::string filetype = bdInfo["filetype"];
  std::string filename = bdInfo["filename"];
  int format = bdInfo["format"];
  int version = bdInfo["version"];
  int numSub = bdInfo["NumSamples"];
  std::vector<int> groupSizes = Rcpp::as<std::vector<int> >(bdInfo["GroupSizes"]);
  Rcpp::DataFrame samples = Rcpp::as<Rcpp::DataFrame>(bdInfo["Samples"]);
  std::vector<std::string> fid = Rcpp::as<std::vector<std::string> >(samples["FID"]);
  std::vector<std::string> sid = Rcpp::as<std::vector<std::string> >(samples["SID"]);
  int numSNPs = bdInfo["NumSNPs"];
  Rcpp::DataFrame snps = Rcpp::as<Rcpp::DataFrame>(bdInfo["SNPs"]);
  std::vector<std::string> snpID = Rcpp::as<std::vector<std::string> >(snps["SNPID"]);
  std::vector<std::string> chromosome = Rcpp::as<std::vector<std::string> >(snps["Chromosome"]);
  std::vector<std::string> ref = Rcpp::as<std::vector<std::string> >(snps["Reference"]);
  std::vector<std::string> alt = Rcpp::as<std::vector<std::string> >(snps["Alternate"]);
  std::vector<int> loc = Rcpp::as<std::vector<int> >(snps["Location"]);
  Rcpp::DataFrame snpInfo = Rcpp::as<Rcpp::DataFrame>(bdInfo["SNPInfo"]);
  std::vector<std::vector<double> > aaf, maf, avgCall, rSq;
  std::vector<double> dblVec;
  std::vector<int>::iterator intIt;
  std::vector<int> index = Rcpp::as< std::vector<int> >(bdInfo["Indices"]);
  std::streampos snpIndex;
  std::vector<std::streampos> posIndex;
  std::vector<std::streampos>::iterator spIt;
  CBDoseMiniReader *miniReader = NULL;
  CBDoseWriter *bdw = NULL;

  if (format == newFormat && newVersion == version) {
    Rcpp::Rcerr << "Old format and new format are the same" << std::endl;
    return retVal;
  }
  if (format < 4)
    miniReader = new CBDoseMiniReader1(filename, numSub, numSNPs);
  else
    miniReader = new CBDoseMiniReader4(filename);

  std::cout << "Created mini readers" << std::endl;
  snpIndex = 0;
  if (index.size() != 0) {
    posIndex.resize(index.size());
    for(intIt = index.begin(), spIt = posIndex.begin(); intIt != index.end(); ++intIt, ++spIt) {
      snpIndex += *intIt;
      *spIt = snpIndex;
      std::cout << '\t' << *spIt;
    }
    std::cout << std::endl;
//    Rcpp::Rcout << "Batching" << std::endl;
    miniReader->ChunkIt(posIndex);
  }

  if (!miniReader->good()) {
    Rcpp::Rcerr << "Error opening binary dosage file" << std::endl;
    delete miniReader;
    return retVal;
  }

  if (newFormat < 4)
    bdw = new CBDoseWriter1(newFile, famFile, mapFile, newFormat, newVersion);
  else
    bdw = new CBDoseWriter4(newFile, newFormat, newVersion);

  if (!bdw->good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }

  if (newFormat == 4 && format == 4) {
    dblVec = Rcpp::as<std::vector<double> >(snpInfo["AAF"]);
    if (dblVec[0] == dblVec[0]) {
      aaf.resize(1);
      aaf[0] = dblVec;
    }
    dblVec = Rcpp::as<std::vector<double> >(snpInfo["MAF"]);
    if (dblVec[0] == dblVec[0]) {
      maf.resize(1);
      maf[0] = dblVec;
    }
    dblVec = Rcpp::as<std::vector<double> >(snpInfo["AvgCall"]);
    if (dblVec[0] == dblVec[0]) {
      avgCall.resize(1);
      avgCall[0] = dblVec;
    }
    dblVec = Rcpp::as<std::vector<double> >(snpInfo["RSquared"]);
    if (dblVec[0] == dblVec[0]) {
      rSq.resize(1);
      rSq[0] = dblVec;
    }
  }

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

    std::cout << "Before writing" << std::endl;
    for (intIt = loc.begin(); intIt != loc.end(); ++intIt) {
      if (intIt == loc.begin())
        miniReader->GetFirst();
      else
        miniReader->GetNext();
      if (bdw->WriteGeneticData(miniReader->Dosage(), miniReader->P0(), miniReader->P1(), miniReader->P2())) {
        Rcpp::Rcerr << "Error genetic data" << std::endl;
        break;
      }
    }
    if (intIt != loc.end())
      break;
  } while (0);
  if (!bdw->good()) {
    Rcpp::Rcerr << "Error writing file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }
  if (bdw->Finalize()) {
    Rcpp::Rcerr << "Error finalizing file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }

  retVal = Rcpp::List::create(Rcpp::Named("SID") = sid,
                              Rcpp::Named("SNPs") = snpID,
                              Rcpp::Named("Dosage") = miniReader->Dosage());
  if (miniReader)
    delete miniReader;
  if (bdw)
    delete bdw;
  return retVal;
}

// [[Rcpp::export]]
int BDConvertVCFC(const Rcpp::List &vcfInfo, const std::string &newFile, const std::string &famFile,
                         const std::string &mapFile, int newFormat, int newVersion) {
  int retVal = 1;
  std::string filetype = vcfInfo["filetype"];
  std::string filename = vcfInfo["filename"];
  int numSub = vcfInfo["NumSamples"];
  Rcpp::DataFrame samples = Rcpp::as<Rcpp::DataFrame>(vcfInfo["Samples"]);
  std::vector<std::string> fid = Rcpp::as<std::vector<std::string> >(samples["FID"]);
  std::vector<std::string> sid = Rcpp::as<std::vector<std::string> >(samples["SID"]);
  int numSNPs = vcfInfo["NumSNPs"];
  Rcpp::DataFrame snps = Rcpp::as<Rcpp::DataFrame>(vcfInfo["SNPs"]);
  std::vector<std::string> snpID = Rcpp::as<std::vector<std::string> >(snps["SNPID"]);
  std::vector<std::string> chromosome = Rcpp::as<std::vector<std::string> >(snps["Chromosome"]);
  std::vector<std::string> ref = Rcpp::as<std::vector<std::string> >(snps["Reference"]);
  std::vector<std::string> alt = Rcpp::as<std::vector<std::string> >(snps["Alternate"]);
  std::vector<int> loc = Rcpp::as<std::vector<int> >(snps["Location"]);
  int startDataI = vcfInfo["StartData"];
  std::streampos startData = startDataI;
  std::vector<std::vector<double> > aaf, maf, avgCall, rSq;
  std::vector<double> dblVec;
  std::vector<int>::iterator intIt;
  std::vector<int> index = Rcpp::as< std::vector<int> >(vcfInfo["Indices"]);
  std::streampos snpIndex;
  std::vector<std::streampos> posIndex;
  std::vector<std::streampos>::iterator spIt;
  CVCFMiniReader *miniReader = NULL;
  CBDoseWriter *bdw = NULL;
  std::vector<int> groupSizes;
  int i;

  miniReader = new CVCFMiniReader(filename, numSub, numSNPs, startData);
  if (!miniReader->good()) {
    Rcpp::Rcerr << "Error opening vcf file" << std::endl;
    delete miniReader;
    return retVal;
  }

  if (newFormat < 4)
    bdw = new CBDoseWriter1(newFile, famFile, mapFile, newFormat, newVersion);
  else
    bdw = new CBDoseWriter4(newFile, newFormat, newVersion);

  snpIndex = 0;
  if (index.size() != 0) {
    posIndex.resize(index.size());
    for(intIt = index.begin(), spIt = posIndex.begin(); intIt != index.end(); ++intIt, ++spIt) {
      snpIndex += *intIt;
      *spIt = snpIndex;
    }
    Rcpp::Rcout << "Batching\t" << posIndex.size() << std::endl;
    miniReader->ChunkIt(posIndex);
    Rcpp::Rcout << "After batching" << std::endl;
  }

  if (!bdw->good()) {
    Rcpp::Rcerr << "Unable to open output file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }
  Rcpp::Rcout << "Before do" << std::endl;
  do {
    groupSizes.push_back(numSub);
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

    i = 1;
    for (intIt = loc.begin(); intIt != loc.end(); ++intIt, ++i) {
      Rcpp::Rcout << i << std::endl;
      if (intIt == loc.begin())
        miniReader->GetFirst();
      else
        miniReader->GetNext();
      if (!miniReader->good()) {
        Rcpp::Rcerr << "Error reading VCF file for SNP:\n" << i << std::endl;
        delete miniReader;
        delete bdw;
        return retVal;
      }
      Rcpp::Rcout << miniReader->Dosage()[0] << '\t' << miniReader->P0()[0] << '\t' << miniReader->P1()[0] << '\t' << miniReader->P2()[0] << std::endl;
      if (bdw->WriteGeneticData(miniReader->Dosage(), miniReader->P0(), miniReader->P1(), miniReader->P2())) {
        Rcpp::Rcout << "Error writing genetic data for SNP:\n" << i << std::endl;
        break;
      }
    }
  } while (0);
  if (!bdw->good()) {
    Rcpp::Rcerr << "Error writing file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }
  if (bdw->Finalize()) {
    Rcpp::Rcerr << "Error finalizing file" << std::endl;
    delete miniReader;
    delete bdw;
    return retVal;
  }

  retVal = 0;
  if (miniReader)
    delete miniReader;
  if (bdw)
    delete bdw;
  return retVal;
}
