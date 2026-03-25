#***************************************************************************#
#                                                                           #
#                      Subsetting Binary Dosage files                       #
#                                                                           #
#***************************************************************************#

#' Subset a binary dosage file
#'
#' Creates a new Format 5 binary dosage file containing a subset of the SNPs
#' and/or subjects from an existing binary dosage file. The input file may be
#' in any format (1-5). At least one filtering criterion must be supplied, and
#' all supplied criteria must be met for a SNP or subject to be retained.
#'
#' @param bdfiles Vector of file names for the input binary dosage file.
#'   Format 4 files require one file name. Formats 1, 2, and 3 require three
#'   file names: the binary dosage file, the family file, and the map file.
#'   Format 5 files require two file names: the .bdose file and the .bdinfo
#'   file.
#' @param bdose_file Path for the output .bdose file.
#' @param bdinfo_file Path for the output .bdinfo file.
#' @param minmaf Minimum minor allele frequency. SNPs whose MAF (computed over
#'   the retained subjects) is below this value are excluded. Must be a single
#'   numeric value between 0 and 0.5.
#' @param locations Integer or numeric vector of SNP base-pair locations to
#'   retain. Cannot be used together with \code{startloc} and \code{endloc}.
#' @param startloc Start of the location range to retain (inclusive). Must be
#'   used together with \code{endloc}. Cannot be used together with
#'   \code{locations}.
#' @param endloc End of the location range to retain (inclusive). Must be used
#'   together with \code{startloc}. Cannot be used together with
#'   \code{locations}.
#' @param subjectids Character vector of subject IDs to retain.
#'
#' @return NULL (invisibly)
#' @export
#'
#' @examples
#' bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo  <- getbdinfo(bdfile)
#' bdose_file  <- tempfile(fileext = ".bdose")
#' bdinfo_file <- tempfile(fileext = ".bdinfo")
#' subsetbd(bdfiles     = bdfile,
#'          bdose_file  = bdose_file,
#'          bdinfo_file = bdinfo_file,
#'          subjectids  = bdinfo$samples$sid[1:30])
subsetbd <- function(bdfiles, bdose_file, bdinfo_file,
                     minmaf = NULL, locations = NULL,
                     startloc = NULL, endloc = NULL,
                     subjectids = NULL) {

  # --- Required argument checks ---
  if (missing(bdfiles))
    stop("No input binary dosage files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (missing(bdinfo_file))
    stop("No output .bdinfo file specified")

  # --- At least one filter required ---
  if (is.null(minmaf) && is.null(locations) &&
      is.null(startloc) && is.null(endloc) && is.null(subjectids))
    stop("At least one of minmaf, locations, startloc/endloc, or subjectids must be specified")

  # --- Location filter mutual exclusivity ---
  if (!is.null(locations) && (!is.null(startloc) || !is.null(endloc)))
    stop("locations and startloc/endloc cannot both be specified")

  # --- startloc/endloc must be paired ---
  if (!is.null(startloc) && is.null(endloc))
    stop("endloc must be specified when startloc is provided")
  if (is.null(startloc) && !is.null(endloc))
    stop("startloc must be specified when endloc is provided")
  if (!is.null(startloc) && startloc > endloc)
    stop("startloc must be less than or equal to endloc")

  # --- minmaf validation ---
  if (!is.null(minmaf)) {
    if (!is.numeric(minmaf) || length(minmaf) != 1L)
      stop("minmaf must be a single numeric value")
    if (minmaf < 0 || minmaf > 0.5)
      stop("minmaf must be between 0 and 0.5")
  }

  # --- Load input bdinfo (formats 1-4 or format 5) ---
  bdinfo <- tryCatch(
    getbdinfo(bdfiles),
    error = function(e) {
      if (length(bdfiles) == 2L) {
        tryCatch(
          getbd5info(bdose_file = bdfiles[1L], bdinfo_file = bdfiles[2L]),
          error = function(e2) stop(conditionMessage(e))
        )
      } else {
        stop(conditionMessage(e))
      }
    }
  )

  n_samp_all <- nrow(bdinfo$samples)
  n_snps_all <- nrow(bdinfo$snps)

  # --- Subject filter ---
  if (!is.null(subjectids)) {
    samp_idx <- which(bdinfo$samples$sid %in% subjectids)
    if (length(samp_idx) == 0L)
      stop("No subjects match the specified subjectids")
    unmatched <- setdiff(subjectids, bdinfo$samples$sid)
    if (length(unmatched) > 0L)
      warning(length(unmatched), " subject ID(s) not found in the file and were ignored",
              call. = FALSE)
  } else {
    samp_idx <- seq_len(n_samp_all)
  }
  n_samp_out <- length(samp_idx)

  # --- SNP location pre-filter (no data read required) ---
  snp_locs <- bdinfo$snps$location
  if (!is.null(locations)) {
    loc_keep <- snp_locs %in% locations
  } else if (!is.null(startloc)) {
    loc_keep <- snp_locs >= startloc & snp_locs <= endloc
  } else {
    loc_keep <- rep(TRUE, n_snps_all)
  }

  # --- Write output .bdose, applying SNP and MAF filters in one pass ---
  bdose_con <- file(bdose_file, open = "wb")
  on.exit(close(bdose_con), add = TRUE)
  writeBin(.bdose5_magic, bdose_con)

  snp_keep   <- logical(n_snps_all)
  indices    <- numeric(n_snps_all)
  current_pos <- 4.0
  n_kept     <- 0L

  for (i in seq_len(n_snps_all)) {
    if (!loc_keep[i]) next

    snp <- getsnp(bdinfo, i, dosageonly = FALSE)

    # MAF filter computed over retained subjects
    if (!is.null(minmaf)) {
      aaf <- mean(snp$dosage[samp_idx], na.rm = TRUE) / 2.0
      maf <- min(aaf, 1.0 - aaf)
      if (is.nan(maf) || maf < minmaf) next
    }

    n_kept       <- n_kept + 1L
    snp_keep[i]  <- TRUE

    ds <- snp$dosage[samp_idx]
    if (all(is.na(snp$p0))) {
      gp <- rep(NA_real_, 3L * n_samp_out)
    } else {
      gp <- as.numeric(rbind(snp$p0[samp_idx], snp$p1[samp_idx], snp$p2[samp_idx]))
    }

    block          <- compress_snp_block(ds, gp)
    indices[n_kept] <- current_pos
    writeBin(block, bdose_con)
    current_pos <- current_pos + length(block)
  }

  if (n_kept == 0L)
    stop("No SNPs pass the specified filters")

  # --- Build and save output bdinfo ---
  out_snps <- bdinfo$snps[snp_keep, , drop = FALSE]
  rownames(out_snps) <- NULL

  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = bdinfo$usesfid,
    samples        = bdinfo$samples[samp_idx, , drop = FALSE],
    onechr         = length(unique(out_snps$chromosome)) == 1L,
    snpidformat    = bdinfo$snpidformat,
    snps           = out_snps,
    snpinfo        = lapply(bdinfo$snpinfo, function(x) x[snp_keep]),
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp_out),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices[seq_len(n_kept)]
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, bdinfo_file)

  invisible(NULL)
}
