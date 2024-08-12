#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <Rcpp.h>
using namespace Rcpp;

//***************************************************************************//
//                           Constants                                       //
//***************************************************************************//

// Magic word for binary dosage files
extern const int MAGICWORD;
// Format ID stored in file
extern const std::vector<std::vector<int> > FORMAT;

// Values that dosages and genetic probabilities are multiplied by to change
// them to short integers

extern const int NUMBEROFBASES;
// 0x7ffe is 32,767 or 2^16 - 1
// 0xfff3 is 65,534 or 2^32 - 1
// 0x2710 is 10,000
extern const unsigned short USBASE[];

// Values the short integers are multiplied by to get dosage and genetic
// probabilities
extern const float FBASE[];

// Tolerance for writing data
extern const double WRITINGTOLERANCE = 1e-6;

// Routine to convert a vector of double values to a vector of unsigned shorts
// Parameter x - is the vector of doubles to convert
// Parameter us - vector of unsigned shorts to store converted values
// Parameter base - Index of USBASE to use as base
// NOTE: missing values are coded as 0xffff
void FloatToUShort(const std::vector<float> &x,
                   std::vector<uint16_t> &us,
                   const int base) {
  uint16_t r1, r2;
  uint16_t *ps1;
  float x1, x2;
  uint32_t i;

  r1 = r2 = 0;
  x1 = x2 = 0.;
  ps1 = (unsigned short *)&us[0];
  for (i = 0; i < x.size(); ++i, ++ps1) {
    if (x[i] != x[i]) {
      // Missing
      *ps1 = 0xffff;
    } else {
      r1 = (unsigned short)(x[i] * USBASE[base]);
      // The following section checks if r1 or (r1 -1) or (r1 + 1)
      // gives the closest value to the double passed
      x1 = r1 * FBASE[base];
      if (r1 * FBASE[base] < x[i])
        r2 = r1 + 1;
      else
        r2 = r1 - 1;
      x2 = r2 * FBASE[base];
      *ps1 = (fabs(x[i] - x1) < fabs(x[i] - x2)) ? r1 : r2;
    }
  }
}

uint16_t FloatToUShort(const float x, const int base) {
  uint16_t r1, r2;
  float x1, x2;

  r1 = r2 = 0;
  x1 = x2 = 0.;

  r1 = (unsigned short)(x * USBASE[base]);
  // The following section checks if r1 or (r1 -1) or (r1 + 1)
  // gives the closest value to the double passed
  x1 = r1 * FBASE[base];
  if (r1 * FBASE[base] < x)
    r2 = r1 + 1;
  else
    r2 = r1 - 1;
  x2 = r2 * FBASE[base];

  return (fabs(x - x1) < fabs(x - x2)) ? r1 : r2;
}


std::vector<float> getFloatsFromTag(bcf_hdr_t *header, bcf1_t *line, const char *tag) {
  int nret = 0;
  float *values = nullptr;
  int n = bcf_get_format_float(header, line, tag, &values, &nret);

  std::vector<float> result;
  if (n > 0) {
    result.assign(values, values + nret);
  } else {
    result.assign(nret, std::numeric_limits<float>::quiet_NaN());
  }

  delete values;
  return result;
}

int getDSFloats(bcf_hdr_t *header, bcf1_t *line, std::vector<float> &d) {
  int nret = 0;
  float *values = nullptr;
  int n = bcf_get_format_float(header, line, "DS", &values, &nret);

  if (n > 0) {
    if (abs(n) != d.size()) {
      Rcpp::Rcerr << "Wrong number of dosages read " << n << '\t' << d.size() << std::endl;
      delete values;
      return 1;
    }
    d.assign(values, values + nret);
  } else {
    d.assign(d.size(), std::numeric_limits<float>::quiet_NaN());
  }

  delete values;
  return 0;
}

int getGPFloats(bcf_hdr_t *header, bcf1_t *line, std::vector<float> &p0,
                std::vector<float> &p1, std::vector<float> &p2) {
  int nret = 0;
  float *pp0, *pp1, *pp2, *pv;
  float *values = nullptr;
  int n = bcf_get_format_float(header, line, "GP", &values, &nret);

  if (n > 0) {
    if (abs(n) != 3*p0.size()) {
      Rcpp::Rcerr << "Wrong number of dosages read " << n << '\t' << 3*p0.size() << std::endl;
      delete values;
      return 1;
    }
    pp0 = &p0[0];
    pp1 = &p1[0];
    pp2 = &p2[0];
    pv = values;
    for (unsigned int i = 0; i < p0.size(); ++i, ++pp0, ++pp1, ++pp2) {
      *pp0 = *pv;
      ++pv;
      *pp1 = *pv;
      ++pv;
      *pp2 = *pv;
      ++pv;
    }
  } else {
    p0.assign(p0.size(), std::numeric_limits<float>::quiet_NaN());
    p1.assign(p1.size(), std::numeric_limits<float>::quiet_NaN());
    p2.assign(p2.size(), std::numeric_limits<float>::quiet_NaN());
  }

  delete values;
  return 0;
}

void CompressGPData(const std::vector<float> &dosage,
                    const std::vector<float> &p0,
                    const std::vector<float> &p1,
                    const std::vector<float> &p2,
                    std::vector<uint16_t> &us,
                    int &datasize) {
  unsigned int i;
  int additionallength;
  uint16_t *usdose, *usadd;

  FloatToUShort(dosage, us, 2);

  additionallength = 0;

  usdose = (uint16_t *)&us[0];
  usadd = usdose + dosage.size();
  for (i = 0; i < dosage.size(); ++i, ++usdose) {
    if (dosage[i] != dosage[i])
      continue;
    if (p0[i] != p0[i] || p1[i] != p1[i] || p2[i] != p2[i]) {
      *usdose |= 0x8000;
      *usadd = 0xffff;
      ++usadd;
      ++additionallength;
    } else if (fabs(p0[i] + p1[i] + p2[i] - 1.) < WRITINGTOLERANCE
                 && fabs(p1[i] + p2[i] + p2[i] - dosage[i]) < WRITINGTOLERANCE) {
      if (p0[i] != 0 && p2[i] != 0) {
        *usdose |= 0x8000;
        *usadd = FloatToUShort(p1[i], 2);
        ++usadd;
        ++additionallength;
      }
    } else {
      *usdose |= 0x8000;
      *usadd = 0x8000 | FloatToUShort(p1[i], 2);
      ++usadd;
      *usadd = FloatToUShort(p0[i], 2);
      ++usadd;
      *usadd = FloatToUShort(p2[i], 2);
      ++usadd;
      additionallength += 3;
    }
  }
  datasize = additionallength + dosage.size();
}

std::vector<Bytef> compressData(const std::vector<uint16_t>& data) {
  uLongf compressedSize = compressBound(data.size() * sizeof(uint16_t));
  std::vector<Bytef> compressedData(compressedSize);

  compress(compressedData.data(), &compressedSize,
           reinterpret_cast<const Bytef*>(data.data()),
           data.size() * sizeof(uint16_t));

  compressedData.resize(compressedSize);
  return compressedData;
}

void writeCompressedData(std::ofstream &outfile, std::vector<Bytef> &compressedData) {
  size_t dataSize = compressedData.size();
  //  outfile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));
  // Write the compressed data
  outfile.write(reinterpret_cast<const char*>(compressedData.data()), dataSize);
}
// [[Rcpp::export]]
Rcpp::List writeVectorsToFile(CharacterVector bd_file_name0,
                              CharacterVector vcf_name0,
                              int dosageOnly)
{
  std::string bd_file_name_str = Rcpp::as<std::string>(bd_file_name0[0]);
  const char* bd_file_name = bd_file_name_str.c_str();
  std::string vcf_name_str = Rcpp::as<std::string>(vcf_name0[0]);
  const char* vcf_name = vcf_name_str.c_str();
  std::vector<std::string> chrom;
  std::vector<int> pos;
  std::vector<std::string> id;
  std::vector<std::string> ref;
  std::vector<std::string> alt;
  const char* c;
  std::string c_s;
  std::string id_s;
  std::string ref_s;
  std::string alt_s;
  std::string alt_s2;
  std::vector<float> gp_values;
  std::vector<float> dosage, p0, p1, p2;
  std::vector<uint16_t> udosage;
  std::vector<uint16_t> udgp;
  std::vector<Bytef> cudgp;
  int usize = 0;
  std::vector<int> dsize;

  htsFile *inf = bcf_open(vcf_name, "r");

  bcf_hdr_t *hdr = bcf_hdr_read(inf);

  bcf1_t *rec = bcf_init();
  int nsamples = bcf_hdr_nsamples(hdr);
  dosage.resize(nsamples, 0);
  p0.resize(nsamples,0);
  p1.resize(nsamples,0);
  p2.resize(nsamples,0);
  udosage.resize(4*nsamples, 0);
  udgp.reserve(4*nsamples);
  size_t dataSize;

  std::ofstream outfile(bd_file_name, std::ios::binary);
  if (!outfile)
  {
    Rcpp::Rcerr << "Error opening file for writing!" << std::endl;
    return NULL;
  }

  Rcpp::CharacterVector sid;
  for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    sid.push_back(hdr->samples[i]);

  outfile.write((char *)&MAGICWORD, sizeof(int));
  if (dosageOnly)
    outfile.write((char *)&FORMAT[4][0], sizeof(int));
  else
    outfile.write((char *)&FORMAT[4][1], sizeof(int));

  usize = nsamples;
  while (bcf_read(inf, hdr, rec) == 0)
  {
    bcf_unpack(rec, BCF_UN_STR);
    bcf_unpack(rec, BCF_UN_INFO);

    //Process CHROM
    c = bcf_hdr_id2name(hdr, rec->rid);
    c_s = c;
    chrom.push_back(c_s);
    //Process POS
    pos.push_back(rec->pos+1);
    //Process ID
    id_s = rec->d.id;
    id.push_back(id_s);
    //process REF
    ref_s = rec->d.allele[0];
    ref.push_back(ref_s);
    //process ALT
    if(rec->n_allele > 1) {
      alt_s = rec->d.allele[1];
    } else {
      alt_s = "None";
    }
    //printf("%lu\n", (unsigned long)rec->n_allele);
    if(rec->n_allele > 2){
      for (int i=2; i<rec->n_allele; ++i){
        alt_s = alt_s + ' ';
        alt_s2 = rec->d.allele[i];
        alt_s = alt_s + alt_s2;
      }
    }
    alt.push_back(alt_s);

    getDSFloats(hdr, rec, dosage);
    FloatToUShort(dosage, udosage, 2);
    if (!dosageOnly) {
      getGPFloats(hdr, rec, p0, p1, p2);
      CompressGPData(dosage, p0, p1, p2, udosage, usize);
    }
    udgp.assign(udosage.begin(), udosage.begin() + usize);
    //    outfile.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));
    cudgp = compressData(udgp);
    writeCompressedData(outfile, cudgp);
    dataSize = cudgp.size();
    dsize.push_back((int)dataSize);

    //    Rcpp::Rcout << usize << '\t' << udgp.size()
    //                << '\t' << cudgp.size() << std::endl;
    //    Rcpp::Rcout << udgp[0] << '\t' << (udgp[0] & 0x7fff) << std::endl;
    //    if(1)
    //      break;
  }

  bcf_destroy(rec);
  bcf_hdr_destroy(hdr);
  if (hts_close(inf))
    Rcpp::Rcerr << "Error closing VCF file" << std::endl;

  outfile.close();
  Rcpp::LogicalVector onechr = true;
  Rcpp::List additionalinfo = Rcpp::List::create(Rcpp::Named("format") = 5,
                                                 Rcpp::Named("subformat") = 2 - dosageOnly,
                                                 Rcpp::Named("headersize") = 8,
                                                 Rcpp::Named("numgroups") = 1,
                                                 Rcpp::Named("groups") = nsamples);
  additionalinfo.attr("class") = "bdose-info";
  Rcpp::List snpinfo = NULL;
  std::vector<double> indices;
  indices.resize(dsize.size());
  indices[0] = 8;
  for (size_t i = 1; i < dsize.size(); ++i)
    indices[i] = indices[i - 1] + dsize[i - 1];
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("filename") = bd_file_name0,
                                         Rcpp::Named("usesfid") = false,
                                         Rcpp::Named("samples") = Rcpp::DataFrame::create(
                                           Rcpp::Named("samples") = sid),
                                           Rcpp::Named("onechr") = onechr,
                                           Rcpp::Named("snpidformat") = 0,
                                           Rcpp::Named("snps") = Rcpp::DataFrame::create(
                                             Rcpp::Named("CHROM") = chrom,
                                             Rcpp::Named("POS") = pos,
                                             Rcpp::Named("ID") = id,
                                             Rcpp::Named("REF") = ref,
                                             Rcpp::Named("ALT") = alt),
                                             Rcpp::Named("snpinfo") = snpinfo,
                                             Rcpp::Named("additionalinfo") = additionalinfo,
                                             Rcpp::Named("datasize") = dsize,
                                             Rcpp::Named("indices") = indices);
  result.attr("class") = "genetic-info";
  return result;
}

