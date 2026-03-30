#***************************************************************************#
#                                                                           #
#               Format 5 Binary Dosage - Write Functions                    #
#                                                                           #
# File pair:                                                                #
#   .bdose  - compressed dosage/genotype probability data                   #
#   .bdi    - SNP metadata, sample IDs, and SNP byte offsets (RDS)          #
#                                                                           #
# Encoding:                                                                 #
#   Values in [0, 2] stored as unsigned short = round(value * 10000)        #
#   Range 0-20000; 0xffff (stored as -1L) = missing                        #
#   DS and GP concatenated per SNP, then gzip-compressed                    #
#                                                                           #
#***************************************************************************#


#***************************************************************************#
#                        Internal constants                                 #
#***************************************************************************#

# Magic bytes at the start of every .bdose Format 5 file:
# byte 0: 0x01 (format major version)
# byte 1: 0x00
# byte 2: 0x05 (format 5)
# byte 3: 0x00
.bdose5_magic <- as.raw(c(0x01, 0x00, 0x05, 0x00))


#***************************************************************************#
#                        Internal helpers                                   #
#***************************************************************************#

# Convert a numeric vector (values in [0, 2], NA = missing) to a raw vector
# of little-endian unsigned 16-bit integers.
# Encoding: stored = round(x * 10000), range 0-20000.
# Missing (NA) is stored as 0xffff, written as -1L in signed int16.
vals_to_ushort_raw <- function(x) {
  stored <- as.integer(round(x * 10000))
  stored[is.na(x)] <- -1L
  writeBin(stored, raw(), size = 2L, endian = "little")
}

# Gzip-compress the DS and GP data for one SNP into a raw vector.
# ds: numeric vector, length n_samples, values in [0, 2]
# gp: numeric vector, length 3 * n_samples, interleaved P(0/0), P(0/1), P(1/1)
compress_snp_block <- function(ds, gp) {
  raw_block <- c(vals_to_ushort_raw(ds), vals_to_ushort_raw(gp))
  memCompress(raw_block, type = "gzip")
}

# Save a Format 5 .bdi file as an RDS object with class "genetic-info".
# The output path is derived as paste0(bdose_file, ".bdi").
#
# Parameters:
#   bdose_file    - path to the associated .bdose file (stored in $filename)
#   samples_sid   - character vector of sample subject IDs
#   snps_df       - data.frame with columns: chromosome, location, snpid,
#                   reference, alternate
#   indices       - numeric vector of byte offsets in the .bdose file
#   snpidformat   - integer; resolved snpidformat value to store in bdinfo
#   snpinfo_list  - optional named list of per-SNP annotation vectors
#                   (e.g. list(aaf = ..., maf = ..., rsq = ...))
write_bdinfo5 <- function(bdose_file, samples_sid,
                          snps_df, indices, snpidformat = 0L,
                          snpinfo_list = list()) {
  bdinfo_file <- paste0(bdose_file, ".bdi")
  n_samp <- length(samples_sid)
  n_snps <- nrow(snps_df)

  samples_df <- data.frame(
    fid = rep("", n_samp),
    sid = samples_sid,
    stringsAsFactors = FALSE
  )

  onechr <- length(unique(snps_df$chromosome)) == 1L

  additional <- list(
    format     = 5,
    subformat  = 1L,
    headersize = 4L,
    numgroups  = 1L,
    groups     = n_samp
  )
  class(additional) <- "bdose-info"

  bdinfo <- list(
    filename       = normalizePath(bdose_file, winslash = "/", mustWork = FALSE),
    usesfid        = FALSE,
    samples        = samples_df,
    onechr         = onechr,
    snpidformat    = snpidformat,
    snps           = snps_df,
    snpinfo        = snpinfo_list,
    additionalinfo = additional,
    datasize       = integer(0L),
    indices        = indices
  )
  class(bdinfo) <- "genetic-info"

  saveRDS(bdinfo, bdinfo_file)
}


#***************************************************************************#
#                        Read functions                                     #
#***************************************************************************#

# Validate a Format 5 .bdose file and load its .bdi metadata.
# Called internally by getbdinfo when a Format 5 magic header is detected.
#
#' @keywords internal
getbd5info <- function(bdose_file) {
  if (missing(bdose_file))
    stop("No .bdose file specified")
  if (!file.exists(bdose_file))
    stop("File not found: ", bdose_file)
  bdinfo_file <- paste0(bdose_file, ".bdi")
  if (!file.exists(bdinfo_file))
    stop("File not found: ", bdinfo_file)

  # Validate the .bdose magic header.
  con <- file(bdose_file, open = "rb")
  magic <- readBin(con, what = raw(), n = 4L)
  close(con)
  if (!identical(magic, .bdose5_magic))
    stop("File does not appear to be a Format 5 .bdose file: ", bdose_file)

  # Load and validate the .bdinfo RDS.
  bdinfo <- readRDS(bdinfo_file)
  if (!is.list(bdinfo))
    stop("Invalid .bdinfo file: not an R list")
  required <- c("snps", "samples", "indices")
  missing_fields <- setdiff(required, names(bdinfo))
  if (length(missing_fields) > 0L)
    stop("Invalid .bdinfo file: missing fields: ",
         paste(missing_fields, collapse = ", "))

  n_snps <- nrow(bdinfo$snps)
  if (length(bdinfo$indices) != n_snps)
    stop("Invalid .bdinfo: indices length (", length(bdinfo$indices),
         ") does not match SNP count (", n_snps, ")")

  # Overwrite filename with the caller-supplied (normalized) path so the
  # object remains valid even if files have been moved.
  bdinfo$filename <- normalizePath(bdose_file, winslash = "/")
  class(bdinfo) <- "genetic-info"
  bdinfo
}

# Decode a raw vector of little-endian signed int16 values back to doubles.
# -1L (stored 0xffff) -> NA; all other values -> x / 10000.
ushort_raw_to_vals <- function(raw_block, n) {
  x <- readBin(raw_block, what = integer(), n = n,
               size = 2L, signed = TRUE, endian = "little")
  result <- as.double(x) / 10000.0
  result[x == -1L] <- NA_real_
  result
}

#' Read a SNP from a Format 5 binary dosage file
#'
#' Seeks to the SNP's compressed block in the .bdose file, decompresses it,
#' and returns the dosage and genotype probabilities for all samples.
#'
#' @param bd5info  Object returned by \code{getbdinfo}.
#' @param snp  The SNP to retrieve: either a 1-based integer index or a
#'   character SNP ID matching a value in \code{bd5info$snps$snpid}.
#'
#' @return A list with four numeric vectors, each of length n_samples:
#' \describe{
#'   \item{dosage}{DS values in \[0, 2\]; NA = missing.}
#'   \item{p0}{P(g=0) values in \[0, 1\]; NA = missing.}
#'   \item{p1}{P(g=1) values in \[0, 1\]; NA = missing.}
#'   \item{p2}{P(g=2) values in \[0, 1\]; NA = missing.}
#' }
#' @export
getbd5snp <- function(bd5info, snp) {
  if (missing(bd5info))
    stop("bd5info missing")
  if (!inherits(bd5info, "genetic-info"))
    stop("bd5info must be an object returned by getbdinfo")
  if (missing(snp))
    stop("No SNP specified")
  if (length(snp) != 1L)
    stop("snp must be of length 1")

  # Resolve SNP to a 1-based integer index.
  if (is.character(snp)) {
    snp_idx <- match(snp, bd5info$snps$snpid)
    if (is.na(snp_idx))
      stop("SNP '", snp, "' not found in bd5info")
  } else {
    snp_idx <- as.integer(floor(snp))
    if (snp_idx < 1L || snp_idx > nrow(bd5info$snps))
      stop("snp index out of range: ", snp_idx)
  }

  n_snps <- nrow(bd5info$snps)
  n_samp <- nrow(bd5info$samples)
  start  <- bd5info$indices[snp_idx]

  # Compressed block size: difference between consecutive offsets, or
  # distance to end-of-file for the last SNP.
  if (snp_idx < n_snps) {
    nbytes <- as.integer(bd5info$indices[snp_idx + 1L] - start)
  } else {
    nbytes <- as.integer(file.info(bd5info$filename)$size - start)
  }

  # Read and decompress the block.
  con <- file(bd5info$filename, open = "rb")
  on.exit(close(con))
  seek(con, where = start, origin = "start")
  compressed <- readBin(con, what = raw(), n = nbytes)
  raw_block  <- memDecompress(compressed, type = "gzip")

  # Decode DS (first n_samp values) and GP (next 3*n_samp values).
  dosage <- ushort_raw_to_vals(raw_block, n_samp)

  gp_raw <- raw_block[seq(2L * n_samp + 1L, length(raw_block))]
  gp     <- ushort_raw_to_vals(gp_raw, 3L * n_samp)

  # Deinterleave GP: layout is [P00_s1, P01_s1, P11_s1, P00_s2, ...]
  idx <- seq_len(n_samp)
  p0  <- gp[3L * (idx - 1L) + 1L]
  p1  <- gp[3L * (idx - 1L) + 2L]
  p2  <- gp[3L * (idx - 1L) + 3L]

  list(dosage = dosage, p0 = p0, p1 = p1, p2 = p2)
}

#' Read a Format 5 SNP into pre-allocated vectors (buffered variant)
#'
#' Like \code{getbd5snp} but writes results into caller-supplied vectors
#' instead of allocating new ones.  Intended for tight loops where thousands
#' of SNPs are read sequentially; pre-allocating the output vectors once
#' avoids repeated memory allocation.
#'
#' The four output vectors \strong{must not} have more than one R binding at
#' the call site (no extra variables pointing to the same object); R's
#' copy-on-modify semantics would otherwise prevent in-place update.
#'
#' @param bd5info  Object returned by \code{getbdinfo}.
#' @param snp  1-based integer index or character SNP ID.
#' @param dosage  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p0  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p1  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p2  Pre-allocated \code{numeric(n_samples)} vector.
#'
#' @return \code{NULL} invisibly.  \code{dosage}, \code{p0}, \code{p1}, and
#'   \code{p2} are updated in place.
#' @export
getbd5snp_buf <- function(bd5info, snp, dosage, p0, p1, p2) {
  if (missing(bd5info))
    stop("bd5info missing")
  if (!inherits(bd5info, "genetic-info"))
    stop("bd5info must be an object returned by getbdinfo")
  if (missing(snp))
    stop("No SNP specified")
  if (length(snp) != 1L)
    stop("snp must be of length 1")

  if (is.character(snp)) {
    snp_idx <- match(snp, bd5info$snps$snpid)
    if (is.na(snp_idx))
      stop("SNP '", snp, "' not found in bd5info")
  } else {
    snp_idx <- as.integer(floor(snp))
    if (snp_idx < 1L || snp_idx > nrow(bd5info$snps))
      stop("snp index out of range: ", snp_idx)
  }

  n_snps <- nrow(bd5info$snps)
  n_samp <- nrow(bd5info$samples)
  start  <- bd5info$indices[snp_idx]

  if (snp_idx < n_snps) {
    nbytes <- as.integer(bd5info$indices[snp_idx + 1L] - start)
  } else {
    nbytes <- as.integer(file.info(bd5info$filename)$size - start)
  }

  con <- file(bd5info$filename, open = "rb")
  on.exit(close(con))
  seek(con, where = start, origin = "start")
  compressed <- readBin(con, what = raw(), n = nbytes)
  raw_block  <- memDecompress(compressed, type = "gzip")

  DecodeFormat5BlockC(raw_block, n_samp, dosage, p0, p1, p2)
  invisible(NULL)
}

#' Open a persistent connection to a Format 5 binary dosage file
#'
#' Opens the .bdose file for reading and returns an object that holds the
#' connection open across multiple calls to \code{getbd5snp_con}.  The
#' connection is closed automatically when the object is garbage-collected
#' or when R exits; call \code{closebd5con} to close it explicitly.
#'
#' @param bd5info  Object returned by \code{getbdinfo}.
#'
#' @return An object of class \code{"bd5con"} to be passed to
#'   \code{getbd5snp_con} and \code{closebd5con}.
#' @export
openbd5con <- function(bd5info) {
  if (missing(bd5info))
    stop("bd5info missing")
  if (!inherits(bd5info, "genetic-info"))
    stop("bd5info must be an object returned by getbdinfo")
  e <- new.env(parent = emptyenv())
  e$xptr <- OpenFormat5FileC(bd5info$filename)
  reg.finalizer(e, function(e) {
    if (!is.null(e$xptr))
      CloseFormat5FileC(e$xptr)
  }, onexit = TRUE)
  class(e) <- "bd5con"
  e
}

#' Close a persistent Format 5 connection
#'
#' Explicitly closes the connection opened by \code{openbd5con}.  Calling
#' this is optional — the finalizer will close it on garbage collection or
#' R exit — but explicit close is preferred to release the file handle
#' promptly.
#'
#' @param bd5con  Object returned by \code{openbd5con}.
#'
#' @return \code{NULL} invisibly.
#' @export
closebd5con <- function(bd5con) {
  if (!inherits(bd5con, "bd5con"))
    stop("bd5con must be an object returned by openbd5con")
  if (!is.null(bd5con$xptr))
    CloseFormat5FileC(bd5con$xptr)
  bd5con$xptr <- NULL
  invisible(NULL)
}

#' Read a Format 5 SNP using a persistent open connection
#'
#' Like \code{getbd5snp_buf} but reuses an already-open file connection
#' instead of opening and closing it on every call.  Use \code{openbd5con}
#' before the loop and \code{closebd5con} (or let the finalizer handle it)
#' after.
#'
#' @param bd5info  Object returned by \code{getbdinfo}.
#' @param snp  1-based integer index or character SNP ID.
#' @param dosage  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p0  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p1  Pre-allocated \code{numeric(n_samples)} vector.
#' @param p2  Pre-allocated \code{numeric(n_samples)} vector.
#' @param bd5con  Object returned by \code{openbd5con}.
#'
#' @return \code{NULL} invisibly.  \code{dosage}, \code{p0}, \code{p1}, and
#'   \code{p2} are updated in place.
#' @export
getbd5snp_con <- function(bd5info, snp, dosage, p0, p1, p2, bd5con) {
  if (missing(bd5info))
    stop("bd5info missing")
  if (!inherits(bd5info, "genetic-info"))
    stop("bd5info must be an object returned by getbdinfo")
  if (missing(snp))
    stop("No SNP specified")
  if (length(snp) != 1L)
    stop("snp must be of length 1")
  if (!inherits(bd5con, "bd5con"))
    stop("bd5con must be an object returned by openbd5con")

  if (is.character(snp)) {
    snp_idx <- match(snp, bd5info$snps$snpid)
    if (is.na(snp_idx))
      stop("SNP '", snp, "' not found in bd5info")
  } else {
    snp_idx <- as.integer(floor(snp))
    if (snp_idx < 1L || snp_idx > nrow(bd5info$snps))
      stop("snp index out of range: ", snp_idx)
  }

  n_snps <- nrow(bd5info$snps)
  n_samp <- nrow(bd5info$samples)
  start  <- bd5info$indices[snp_idx]

  if (snp_idx < n_snps) {
    nbytes <- as.integer(bd5info$indices[snp_idx + 1L] - start)
  } else {
    nbytes <- as.integer(file.info(bd5info$filename)$size - start)
  }

  ReadFormat5SNPC(bd5con$xptr, start, nbytes, n_samp, dosage, p0, p1, p2)
  invisible(NULL)
}

#***************************************************************************#
#                        Main conversion function                           #
#***************************************************************************#

#' Convert a VCF file to Format 5 binary dosage files
#'
#' Reads the DS (dosage) and GP (genotype probabilities) FORMAT fields from a
#' bgzipped, tabix-indexed VCF file — as produced by imputation servers such
#' as the Michigan Imputation Server — and writes a pair of Format 5
#' BinaryDosage files.
#'
#' The .bdose file begins with a 4-byte magic number followed by one
#' gzip-compressed block per SNP.  Each block contains the DS values for all
#' samples followed by the GP values, encoded as unsigned 16-bit integers
#' (round(value * 10000); 0xffff = missing).
#'
#' The .bdi file is an RDS-serialised R list of class \code{"genetic-info"}
#' with the following elements:
#' \describe{
#'   \item{filename}{Path to the associated .bdose file.}
#'   \item{usesfid}{Logical; always FALSE for VCF-sourced files.}
#'   \item{samples}{data.frame with columns \code{fid} (empty) and
#'     \code{sid} (sample IDs).}
#'   \item{onechr}{Logical; TRUE if all SNPs are on a single chromosome.}
#'   \item{snpidformat}{Numeric; resolved SNP ID format (see \code{snpidformat}
#'     parameter).}
#'   \item{snps}{data.frame with columns chromosome, location, snpid,
#'     reference, alternate.}
#'   \item{snpinfo}{Named list of per-SNP annotations requested via
#'     \code{bdoptions}.  Each element is a numeric vector of length equal to
#'     the number of SNPs.  Values are read from the VCF INFO column when
#'     available for the first SNP (AF for aaf, MAF for maf, R2 for rsq);
#'     otherwise they are calculated from the dosage values.}
#'   \item{additionalinfo}{List of class \code{"bdose-info"} with format,
#'     subformat, headersize, numgroups, and groups.}
#'   \item{datasize}{Integer vector of length 0 (unused in Format 5).}
#'   \item{indices}{Numeric vector of byte offsets into .bdose, one per SNP.}
#' }
#'
#' @param vcffile  Path to the bgzipped, tabix-indexed VCF file.
#' @param bdose_file  Path for the output .bdose file. The companion .bdi
#'   metadata file is written to \code{paste0(bdose_file, ".bdi")}.
#' @param region  Optional genomic region string in bcftools format
#'   (e.g. \code{"chr21"} or \code{"chr21:1-5000000"}).  Requires a
#'   tabix index.  Default \code{NULL} processes the entire file.
#' @param snpidformat  Integer controlling how SNP IDs are stored.
#'   \describe{
#'     \item{-1}{Generate IDs as \code{chr:pos:ref:alt}; equivalent to 2 for
#'       Format 5.}
#'     \item{0}{Use the IDs as they appear in the VCF file (default).
#'       Auto-detects format 1 or 2 if all IDs match.}
#'     \item{1}{Store IDs as \code{chr:pos}.  An error is raised if the VCF
#'       already uses \code{chr:pos:ref:alt} format, as information would be
#'       lost.}
#'     \item{2}{Store IDs as \code{chr:pos:ref:alt}.}
#'     \item{3}{Store IDs as \code{chr:pos_ref_alt}.}
#'   }
#'
#' @param bdoptions  Character vector specifying which per-SNP statistics to
#'   store.  Any combination of \code{"aaf"} (alternate allele frequency),
#'   \code{"maf"} (minor allele frequency), and \code{"rsq"} (imputation
#'   r-squared).  For each statistic, the corresponding VCF INFO field is used
#'   when present for the first SNP (AF, MAF, R2 respectively); otherwise the
#'   value is calculated from the dosage data.  Default \code{character(0)}
#'   stores no statistics.
#'
#' @return NULL (invisibly)
#' @importFrom stats var
#' @importFrom vcfppR vcfreader
#' @export
vcftobd <- function(vcffile, bdose_file, region = NULL,
                     snpidformat = 0L, bdoptions = character(0)) {

  if (!file.exists(vcffile))
    stop("VCF file not found: ", vcffile)

  if (!is.null(region) && (!is.character(region) || length(region) != 1L))
    stop("region must be a single character string or NULL")

  if (!is.numeric(snpidformat) || length(snpidformat) != 1L ||
      floor(snpidformat) != snpidformat)
    stop("snpidformat must be an integer value")
  snpidformat <- as.integer(snpidformat)
  if (snpidformat < -1L || snpidformat > 3L)
    stop("snpidformat must be -1, 0, 1, 2, or 3")

  if (!is.character(bdoptions))
    stop("bdoptions must be a character vector")
  invalid_opts <- setdiff(bdoptions, c("aaf", "maf", "rsq"))
  if (length(invalid_opts) > 0L)
    stop("Invalid bdoptions value(s): ",
         paste(invalid_opts, collapse = ", "),
         ". Valid values are: aaf, maf, rsq")

  want_aaf <- "aaf" %in% bdoptions
  want_maf <- "maf" %in% bdoptions
  want_rsq <- "rsq" %in% bdoptions

  # Pre-allocate metadata vectors; will trim to actual SNP count at the end.
  # 2 million SNPs is sufficient for any single chromosome.
  MAX_SNPS    <- 2000000L
  snp_chr     <- character(MAX_SNPS)
  snp_pos     <- integer(MAX_SNPS)
  snp_id      <- character(MAX_SNPS)
  snp_ref     <- character(MAX_SNPS)
  snp_alt     <- character(MAX_SNPS)
  snp_indices <- numeric(MAX_SNPS)
  snp_count   <- 0L

  # Pre-allocate snpinfo accumulation vectors and source flags.
  # Flags (use_*_from_vcf) are set on the first SNP based on INFO availability.
  if (want_aaf) { aaf_vec <- numeric(MAX_SNPS); use_aaf_from_vcf <- FALSE }
  if (want_maf) { maf_vec <- numeric(MAX_SNPS); use_maf_from_vcf <- FALSE }
  if (want_rsq) { rsq_vec <- numeric(MAX_SNPS); use_rsq_from_vcf <- FALSE }

  # current_pos tracks the write position in the .bdose file (numeric to
  # handle files larger than 2 GB without integer overflow).
  current_pos <- 4  # bytes consumed by the 4-byte header

  # Open .bdose and write the magic header.
  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  # Open the VCF reader and retrieve sample IDs.
  if (!is.null(region))
    reader <- vcfppR::vcfreader$new(vcffile, region)
  else
    reader <- vcfppR::vcfreader$new(vcffile)

  samples_sid <- reader$samples()

  # Stream through every variant in the VCF.
  while (reader$variant()) {
    snp_count <- snp_count + 1L

    if (snp_count > MAX_SNPS)
      stop("SNP count exceeds MAX_SNPS (", MAX_SNPS, "). ",
           "Increase MAX_SNPS in vcftobd().")

    # Validate required FORMAT tags and detect INFO field availability on
    # the first SNP.
    if (snp_count == 1L) {
      if (reader$getFormatType("DS") == 0L)
        stop("FORMAT field 'DS' not found in VCF file: ", vcffile)
      if (reader$getFormatType("GP") == 0L)
        stop("FORMAT field 'GP' not found in VCF file: ", vcffile)
      if (want_aaf) {
        val <- tryCatch(reader$infoFloat("AF"),  error = function(e) NA_real_)
        use_aaf_from_vcf <- length(val) > 0L && !is.na(val[1L])
      }
      if (want_maf) {
        val <- tryCatch(reader$infoFloat("MAF"), error = function(e) NA_real_)
        use_maf_from_vcf <- length(val) > 0L && !is.na(val[1L])
      }
      if (want_rsq) {
        val <- tryCatch(reader$infoFloat("R2"),  error = function(e) NA_real_)
        use_rsq_from_vcf <- length(val) > 0L && !is.na(val[1L])
      }
    }

    # Collect SNP metadata.
    snp_chr[snp_count] <- reader$chr()
    snp_pos[snp_count] <- reader$pos()
    snp_id[snp_count]  <- reader$id()
    snp_ref[snp_count] <- reader$ref()
    snp_alt[snp_count] <- reader$alt()

    # Read dosage and genotype probability data.
    ds <- reader$formatFloat("DS")
    gp <- reader$formatFloat("GP")

    # Accumulate snpinfo values.
    if (want_aaf) {
      aaf_vec[snp_count] <- if (use_aaf_from_vcf)
        tryCatch(reader$infoFloat("AF")[1L],  error = function(e) NA_real_)
      else
        mean(ds, na.rm = TRUE) / 2
    }
    if (want_maf) {
      maf_vec[snp_count] <- if (use_maf_from_vcf)
        tryCatch(reader$infoFloat("MAF")[1L], error = function(e) NA_real_)
      else {
        aaf_i <- if (want_aaf) aaf_vec[snp_count] else mean(ds, na.rm = TRUE) / 2
        min(aaf_i, 1 - aaf_i)
      }
    }
    if (want_rsq) {
      rsq_vec[snp_count] <- if (use_rsq_from_vcf)
        tryCatch(reader$infoFloat("R2")[1L],  error = function(e) NA_real_)
      else {
        aaf_i <- if (want_aaf) aaf_vec[snp_count] else mean(ds, na.rm = TRUE) / 2
        denom <- 2 * aaf_i * (1 - aaf_i)
        if (denom > 0) var(ds, na.rm = TRUE) / denom else NA_real_
      }
    }

    # Record the byte offset for this SNP, compress, and write.
    snp_indices[snp_count] <- current_pos
    compressed <- compress_snp_block(ds, gp)
    writeBin(compressed, bdose_con)
    current_pos <- current_pos + length(compressed)
  }

  if (snp_count == 0L)
    stop("No variants found in VCF file: ", vcffile)

  # Trim metadata vectors to actual SNP count.
  chr_vec <- snp_chr[1:snp_count]
  pos_vec <- snp_pos[1:snp_count]
  ref_vec <- snp_ref[1:snp_count]
  alt_vec <- snp_alt[1:snp_count]
  id_vec  <- snp_id[1:snp_count]

  # Apply snpidformat: resolve IDs and the stored format value.
  chrlocid          <- paste(chr_vec, pos_vec, sep = ":")
  chrlocrefaltid    <- paste(chr_vec, pos_vec, ref_vec, alt_vec, sep = ":")
  chrlocrefaltid_us <- paste(chrlocid, ref_vec, alt_vec, sep = "_")

  if (snpidformat == 0L) {
    # Auto-detect: check whether all VCF IDs already match a known format.
    if (all(id_vec == chrlocid))
      snpidformat <- 1L
    else if (all(id_vec == chrlocrefaltid))
      snpidformat <- 2L
    else if (all(id_vec == chrlocrefaltid_us))
      snpidformat <- 3L
    # else stays 0 — IDs are used verbatim
  } else if (snpidformat == 1L) {
    if (all(id_vec == chrlocrefaltid) || all(id_vec == chrlocrefaltid_us))
      stop("snpidformat 1 specified but VCF file uses a format that encodes ",
           "ref/alt; information would be lost")
    id_vec <- chrlocid
  } else if (snpidformat == 2L) {
    id_vec <- chrlocrefaltid
  } else if (snpidformat == 3L) {
    id_vec <- chrlocrefaltid_us
  } else {
    # snpidformat == -1: generate chr:pos:ref:alt, store as format 2.
    id_vec      <- chrlocrefaltid
    snpidformat <- 2L
  }

  # Build the SNP data frame.
  snps_df <- data.frame(
    chromosome = chr_vec,
    location   = pos_vec,
    snpid      = id_vec,
    reference  = ref_vec,
    alternate  = alt_vec,
    stringsAsFactors = FALSE
  )

  # Build the snpinfo list (trim vectors to actual SNP count).
  snpinfo_list <- list()
  if (want_aaf) snpinfo_list$aaf <- aaf_vec[1:snp_count]
  if (want_maf) snpinfo_list$maf <- maf_vec[1:snp_count]
  if (want_rsq) snpinfo_list$rsq <- rsq_vec[1:snp_count]

  # Save the .bdi file.
  write_bdinfo5(bdose_file   = bdose_file,
                samples_sid  = samples_sid,
                snps_df      = snps_df,
                indices      = snp_indices[1:snp_count],
                snpidformat  = snpidformat,
                snpinfo_list = snpinfo_list)

  invisible(NULL)
}
