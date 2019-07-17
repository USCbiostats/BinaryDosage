// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ReadBinaryDosageBaseHeader
Rcpp::List ReadBinaryDosageBaseHeader(std::string& filename);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageBaseHeader(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageBaseHeader(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageHeader3A
Rcpp::List ReadBinaryDosageHeader3A(std::string& filename);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageHeader3A(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageHeader3A(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageHeader3B
Rcpp::List ReadBinaryDosageHeader3B(std::string& filename);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageHeader3B(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageHeader3B(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageHeader4A
Rcpp::List ReadBinaryDosageHeader4A(std::string& filename);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageHeader4A(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageHeader4A(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageHeader4B
Rcpp::List ReadBinaryDosageHeader4B(std::string& filename);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageHeader4B(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageHeader4B(filename));
    return rcpp_result_gen;
END_RCPP
}
// ReadBDIndicesS4
Rcpp::IntegerVector ReadBDIndicesS4(std::string filename, int numSNPs, int indexStart);
RcppExport SEXP _BinaryDosage_ReadBDIndicesS4(SEXP filenameSEXP, SEXP numSNPsSEXP, SEXP indexStartSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    Rcpp::traits::input_parameter< int >::type indexStart(indexStartSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBDIndicesS4(filename, numSNPs, indexStart));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageDataC
int ReadBinaryDosageDataC(std::string& filename, int headersize, int snp, Rcpp::NumericVector& dosage, Rcpp::IntegerVector& usdosage, int base);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageDataC(SEXP filenameSEXP, SEXP headersizeSEXP, SEXP snpSEXP, SEXP dosageSEXP, SEXP usdosageSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type headersize(headersizeSEXP);
    Rcpp::traits::input_parameter< int >::type snp(snpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type usdosage(usdosageSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageDataC(filename, headersize, snp, dosage, usdosage, base));
    return rcpp_result_gen;
END_RCPP
}
// ReadBinaryDosageDataP1P2
int ReadBinaryDosageDataP1P2(std::string& filename, int headersize, int snp, Rcpp::NumericVector& dosage, Rcpp::NumericVector& p0, Rcpp::NumericVector& p1, Rcpp::NumericVector& p2, Rcpp::IntegerVector& us, int base);
RcppExport SEXP _BinaryDosage_ReadBinaryDosageDataP1P2(SEXP filenameSEXP, SEXP headersizeSEXP, SEXP snpSEXP, SEXP dosageSEXP, SEXP p0SEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP usSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type headersize(headersizeSEXP);
    Rcpp::traits::input_parameter< int >::type snp(snpSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type us(usSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadBinaryDosageDataP1P2(filename, headersize, snp, dosage, p0, p1, p2, us, base));
    return rcpp_result_gen;
END_RCPP
}
// OpenVCFFile
Rcpp::List OpenVCFFile(Rcpp::StringVector& rFilename);
RcppExport SEXP _BinaryDosage_OpenVCFFile(SEXP rFilenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type rFilename(rFilenameSEXP);
    rcpp_result_gen = Rcpp::wrap(OpenVCFFile(rFilename));
    return rcpp_result_gen;
END_RCPP
}
// ReadVCFValues
int ReadVCFValues(std::string& filename, Rcpp::NumericMatrix& values, Rcpp::IntegerVector& snps, Rcpp::IntegerVector& subjects, Rcpp::NumericVector& indices, Rcpp::StringVector& formats);
RcppExport SEXP _BinaryDosage_ReadVCFValues(SEXP filenameSEXP, SEXP valuesSEXP, SEXP snpsSEXP, SEXP subjectsSEXP, SEXP indicesSEXP, SEXP formatsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snps(snpsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type subjects(subjectsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type formats(formatsSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadVCFValues(filename, values, snps, subjects, indices, formats));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageBaseHeader
int WriteBinaryDosageBaseHeader(std::string& filename, int format, int subformat);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageBaseHeader(SEXP filenameSEXP, SEXP formatSEXP, SEXP subformatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type format(formatSEXP);
    Rcpp::traits::input_parameter< int >::type subformat(subformatSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageBaseHeader(filename, format, subformat));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageHeader3A
int WriteBinaryDosageHeader3A(std::string& filename, int numSubjects);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageHeader3A(SEXP filenameSEXP, SEXP numSubjectsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSubjects(numSubjectsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageHeader3A(filename, numSubjects));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageHeader3B
int WriteBinaryDosageHeader3B(std::string& filename, std::string& md5samples, std::string& md5SNPs);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageHeader3B(SEXP filenameSEXP, SEXP md5samplesSEXP, SEXP md5SNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string& >::type md5samples(md5samplesSEXP);
    Rcpp::traits::input_parameter< std::string& >::type md5SNPs(md5SNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageHeader3B(filename, md5samples, md5SNPs));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageHeader4A
int WriteBinaryDosageHeader4A(std::string& filename, int numSubjects, int numSNPs);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageHeader4A(SEXP filenameSEXP, SEXP numSubjectsSEXP, SEXP numSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSubjects(numSubjectsSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageHeader4A(filename, numSubjects, numSNPs));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageHeader4B
int WriteBinaryDosageHeader4B(std::string& filename, int numSubjects, int numSNPs);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageHeader4B(SEXP filenameSEXP, SEXP numSubjectsSEXP, SEXP numSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSubjects(numSubjectsSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageHeader4B(filename, numSubjects, numSNPs));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDGroups
int WriteBDGroups(std::string& filename, Rcpp::IntegerVector& groups);
RcppExport SEXP _BinaryDosage_WriteBDGroups(SEXP filenameSEXP, SEXP groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type groups(groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDGroups(filename, groups));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDGroups2
int WriteBDGroups2(std::string& filename, Rcpp::IntegerVector& groups);
RcppExport SEXP _BinaryDosage_WriteBDGroups2(SEXP filenameSEXP, SEXP groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type groups(groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDGroups2(filename, groups));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDFamilyInfoC
int WriteBDFamilyInfoC(std::string& filename, int numSub, std::string& sid, std::string& fid, int numSubLoc, int suboffsetLoc, int snpoffsetLoc);
RcppExport SEXP _BinaryDosage_WriteBDFamilyInfoC(SEXP filenameSEXP, SEXP numSubSEXP, SEXP sidSEXP, SEXP fidSEXP, SEXP numSubLocSEXP, SEXP suboffsetLocSEXP, SEXP snpoffsetLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSub(numSubSEXP);
    Rcpp::traits::input_parameter< std::string& >::type sid(sidSEXP);
    Rcpp::traits::input_parameter< std::string& >::type fid(fidSEXP);
    Rcpp::traits::input_parameter< int >::type numSubLoc(numSubLocSEXP);
    Rcpp::traits::input_parameter< int >::type suboffsetLoc(suboffsetLocSEXP);
    Rcpp::traits::input_parameter< int >::type snpoffsetLoc(snpoffsetLocSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDFamilyInfoC(filename, numSub, sid, fid, numSubLoc, suboffsetLoc, snpoffsetLoc));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDSNPInfoC
int WriteBDSNPInfoC(std::string& filename, int numSNPs, std::string& snpid, std::string& chromosome, Rcpp::IntegerVector& location, std::string& reference, std::string& alternate, Rcpp::NumericVector& aaf, Rcpp::NumericVector& maf, Rcpp::NumericVector& avgCall, Rcpp::NumericVector& rsq, int numSNPloc, int snpOptionsLoc, int snpOffsetLoc, int nextOffsetLoc);
RcppExport SEXP _BinaryDosage_WriteBDSNPInfoC(SEXP filenameSEXP, SEXP numSNPsSEXP, SEXP snpidSEXP, SEXP chromosomeSEXP, SEXP locationSEXP, SEXP referenceSEXP, SEXP alternateSEXP, SEXP aafSEXP, SEXP mafSEXP, SEXP avgCallSEXP, SEXP rsqSEXP, SEXP numSNPlocSEXP, SEXP snpOptionsLocSEXP, SEXP snpOffsetLocSEXP, SEXP nextOffsetLocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    Rcpp::traits::input_parameter< std::string& >::type snpid(snpidSEXP);
    Rcpp::traits::input_parameter< std::string& >::type chromosome(chromosomeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type location(locationSEXP);
    Rcpp::traits::input_parameter< std::string& >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< std::string& >::type alternate(alternateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type aaf(aafSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type avgCall(avgCallSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type rsq(rsqSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPloc(numSNPlocSEXP);
    Rcpp::traits::input_parameter< int >::type snpOptionsLoc(snpOptionsLocSEXP);
    Rcpp::traits::input_parameter< int >::type snpOffsetLoc(snpOffsetLocSEXP);
    Rcpp::traits::input_parameter< int >::type nextOffsetLoc(nextOffsetLocSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDSNPInfoC(filename, numSNPs, snpid, chromosome, location, reference, alternate, aaf, maf, avgCall, rsq, numSNPloc, snpOptionsLoc, snpOffsetLoc, nextOffsetLoc));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDIndexArray3_4
int WriteBDIndexArray3_4(std::string& filename, int numSNPs);
RcppExport SEXP _BinaryDosage_WriteBDIndexArray3_4(SEXP filenameSEXP, SEXP numSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDIndexArray3_4(filename, numSNPs));
    return rcpp_result_gen;
END_RCPP
}
// WriteBDIndexArray4
int WriteBDIndexArray4(std::string& filename, int numSNPs, int indexoffsetLoc, int dosageoffsetloc);
RcppExport SEXP _BinaryDosage_WriteBDIndexArray4(SEXP filenameSEXP, SEXP numSNPsSEXP, SEXP indexoffsetLocSEXP, SEXP dosageoffsetlocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< int >::type numSNPs(numSNPsSEXP);
    Rcpp::traits::input_parameter< int >::type indexoffsetLoc(indexoffsetLocSEXP);
    Rcpp::traits::input_parameter< int >::type dosageoffsetloc(dosageoffsetlocSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBDIndexArray4(filename, numSNPs, indexoffsetLoc, dosageoffsetloc));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryDosageDataC
int WriteBinaryDosageDataC(const std::string& filename, Rcpp::NumericVector& dosage, Rcpp::IntegerVector& us, int base);
RcppExport SEXP _BinaryDosage_WriteBinaryDosageDataC(SEXP filenameSEXP, SEXP dosageSEXP, SEXP usSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type dosage(dosageSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type us(usSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryDosageDataC(filename, dosage, us, base));
    return rcpp_result_gen;
END_RCPP
}
// WriteBinaryP1P2Data
int WriteBinaryP1P2Data(const std::string& filename, Rcpp::NumericVector& p1, Rcpp::NumericVector& p2, Rcpp::IntegerVector& us, int base);
RcppExport SEXP _BinaryDosage_WriteBinaryP1P2Data(SEXP filenameSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP usSEXP, SEXP baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type us(usSEXP);
    Rcpp::traits::input_parameter< int >::type base(baseSEXP);
    rcpp_result_gen = Rcpp::wrap(WriteBinaryP1P2Data(filename, p1, p2, us, base));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BinaryDosage_ReadBinaryDosageBaseHeader", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageBaseHeader, 1},
    {"_BinaryDosage_ReadBinaryDosageHeader3A", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageHeader3A, 1},
    {"_BinaryDosage_ReadBinaryDosageHeader3B", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageHeader3B, 1},
    {"_BinaryDosage_ReadBinaryDosageHeader4A", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageHeader4A, 1},
    {"_BinaryDosage_ReadBinaryDosageHeader4B", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageHeader4B, 1},
    {"_BinaryDosage_ReadBDIndicesS4", (DL_FUNC) &_BinaryDosage_ReadBDIndicesS4, 3},
    {"_BinaryDosage_ReadBinaryDosageDataC", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageDataC, 6},
    {"_BinaryDosage_ReadBinaryDosageDataP1P2", (DL_FUNC) &_BinaryDosage_ReadBinaryDosageDataP1P2, 9},
    {"_BinaryDosage_OpenVCFFile", (DL_FUNC) &_BinaryDosage_OpenVCFFile, 1},
    {"_BinaryDosage_ReadVCFValues", (DL_FUNC) &_BinaryDosage_ReadVCFValues, 6},
    {"_BinaryDosage_WriteBinaryDosageBaseHeader", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageBaseHeader, 3},
    {"_BinaryDosage_WriteBinaryDosageHeader3A", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageHeader3A, 2},
    {"_BinaryDosage_WriteBinaryDosageHeader3B", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageHeader3B, 3},
    {"_BinaryDosage_WriteBinaryDosageHeader4A", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageHeader4A, 3},
    {"_BinaryDosage_WriteBinaryDosageHeader4B", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageHeader4B, 3},
    {"_BinaryDosage_WriteBDGroups", (DL_FUNC) &_BinaryDosage_WriteBDGroups, 2},
    {"_BinaryDosage_WriteBDGroups2", (DL_FUNC) &_BinaryDosage_WriteBDGroups2, 2},
    {"_BinaryDosage_WriteBDFamilyInfoC", (DL_FUNC) &_BinaryDosage_WriteBDFamilyInfoC, 7},
    {"_BinaryDosage_WriteBDSNPInfoC", (DL_FUNC) &_BinaryDosage_WriteBDSNPInfoC, 15},
    {"_BinaryDosage_WriteBDIndexArray3_4", (DL_FUNC) &_BinaryDosage_WriteBDIndexArray3_4, 2},
    {"_BinaryDosage_WriteBDIndexArray4", (DL_FUNC) &_BinaryDosage_WriteBDIndexArray4, 4},
    {"_BinaryDosage_WriteBinaryDosageDataC", (DL_FUNC) &_BinaryDosage_WriteBinaryDosageDataC, 4},
    {"_BinaryDosage_WriteBinaryP1P2Data", (DL_FUNC) &_BinaryDosage_WriteBinaryP1P2Data, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_BinaryDosage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
