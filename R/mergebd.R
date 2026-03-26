#***************************************************************************#
#                                                                           #
#                      Merging Format 5 Binary Dosage files                 #
#                                                                           #
#***************************************************************************#

#' Merge Format 5 binary dosage files
#'
#' Merges two or more Format 5 binary dosage files into a single Format 5
#' output file. The merge type is determined automatically:
#'
#' \itemize{
#'   \item If subject IDs do not overlap across files, a \strong{subject merge}
#'         is performed: the output contains all subjects from every file and
#'         the SNPs common to all files.
#'   \item If SNP IDs do not overlap across files, a \strong{SNP merge} is
#'         performed: the output contains all SNPs from every file and the
#'         subjects common to all files.
#' }
#'
#' If both subject IDs and SNP IDs overlap across files the merge cannot be
#' performed and an error is returned.
#'
#' SNPs are identified by chromosome, position, reference allele, and alternate
#' allele.
#'
#' @param bdose_files Character vector of paths to the input .bdose files.
#'   Must contain at least two entries. The companion .bdi file for each is
#'   expected at \code{paste0(bdose_files[i], ".bdi")}.
#' @param bdose_file Path for the output .bdose file. The companion .bdi
#'   metadata file is written to \code{paste0(bdose_file, ".bdi")}.
#'
#' @return NULL (invisibly)
#' @export
#'
#' @examples
#' bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo_src <- getbdinfo(bdfile)
#'
#' # Create two format 5 files with non-overlapping subjects
#' bdose_a   <- tempfile(fileext = ".bdose")
#' bdose_b   <- tempfile(fileext = ".bdose")
#' bdose_tmp <- tempfile(fileext = ".bdose")
#' updatebd(bdfiles = bdfile, bdose_file = bdose_tmp)
#' subsetbd(bdfiles    = bdose_tmp,
#'          bdose_file = bdose_a,
#'          subjectids = bdinfo_src$samples$sid[1:30])
#' subsetbd(bdfiles    = bdose_tmp,
#'          bdose_file = bdose_b,
#'          subjectids = bdinfo_src$samples$sid[31:60])
#'
#' bdose_out <- tempfile(fileext = ".bdose")
#' mergebd(bdose_files = c(bdose_a, bdose_b),
#'         bdose_file  = bdose_out)
mergebd <- function(bdose_files, bdose_file) {

  # --- Argument validation ---
  if (missing(bdose_files))
    stop("No input .bdose files specified")
  if (missing(bdose_file))
    stop("No output .bdose file specified")
  if (!is.character(bdose_files))
    stop("bdose_files must be a character vector")
  if (length(bdose_files) < 2L)
    stop("At least two input files are required for merging")

  n_files <- length(bdose_files)

  # --- Load all bdinfo objects (getbd5info verifies Format 5) ---
  bdinfos <- vector("list", n_files)
  for (k in seq_len(n_files)) {
    bdinfos[[k]] <- getbd5info(bdose_file = bdose_files[k])
  }

  n_snps <- vapply(bdinfos, function(b) nrow(b$snps),    integer(1L))
  n_samp <- vapply(bdinfos, function(b) nrow(b$samples), integer(1L))

  # SNP key: chr:pos:ref:alt (unambiguous regardless of snpidformat)
  snp_keys <- lapply(bdinfos, function(b)
    paste(b$snps$chromosome, b$snps$location,
          b$snps$reference,  b$snps$alternate, sep = ":"))
  subj_ids <- lapply(bdinfos, function(b) b$samples$sid)

  # --- Determine merge type ---
  subjects_intersect <- anyDuplicated(unlist(subj_ids))  > 0L
  snps_intersect     <- anyDuplicated(unlist(snp_keys)) > 0L

  if (subjects_intersect && snps_intersect)
    stop("Cannot merge: subject IDs and SNP IDs both overlap across files")

  if (!subjects_intersect) {
    # ------------------------------------------------------------------ #
    #  Subject merge: unique subjects across files, SNPs must overlap     #
    # ------------------------------------------------------------------ #
    common_keys <- snp_keys[[1L]]
    for (k in seq(2L, n_files))
      common_keys <- intersect(common_keys, snp_keys[[k]])

    if (length(common_keys) == 0L)
      stop("Cannot merge: no SNPs in common across files for subject merge")

    # Index of each common SNP in each file
    snp_idx <- lapply(snp_keys, function(kv) match(common_keys, kv))

    n_snps_out <- length(common_keys)
    n_samp_out <- sum(n_samp)

    bdose_con <- file(bdose_file, open = "wb")
    on.exit(close(bdose_con), add = TRUE)
    writeBin(.bdose5_magic, bdose_con)

    indices     <- numeric(n_snps_out)
    current_pos <- 4.0

    for (i in seq_len(n_snps_out)) {
      ds_all <- numeric(0)
      gp_all <- numeric(0)

      for (k in seq_len(n_files)) {
        snp <- getsnp(bdinfos[[k]], snp_idx[[k]][i], dosageonly = FALSE)
        ds_all <- c(ds_all, snp$dosage)
        if (all(is.na(snp$p0))) {
          gp_all <- c(gp_all, rep(NA_real_, 3L * n_samp[k]))
        } else {
          gp_all <- c(gp_all,
                      as.numeric(rbind(snp$p0, snp$p1, snp$p2)))
        }
      }

      block       <- compress_snp_block(ds_all, gp_all)
      indices[i]  <- current_pos
      writeBin(block, bdose_con)
      current_pos <- current_pos + length(block)
    }

    # Output SNP info from file 1 (in common-SNP order)
    out_snps <- bdinfos[[1L]]$snps[snp_idx[[1L]], , drop = FALSE]
    rownames(out_snps) <- NULL

    # Output samples: all subjects in file order
    out_samples <- do.call(rbind, lapply(bdinfos, function(b) b$samples))
    rownames(out_samples) <- NULL

  } else {
    # ------------------------------------------------------------------ #
    #  SNP merge: unique SNPs across files, subjects must overlap         #
    # ------------------------------------------------------------------ #
    common_sids <- subj_ids[[1L]]
    for (k in seq(2L, n_files))
      common_sids <- intersect(common_sids, subj_ids[[k]])

    if (length(common_sids) == 0L)
      stop("Cannot merge: no subjects in common across files for SNP merge")

    # Index of each common subject in each file
    samp_idx <- lapply(subj_ids, function(sv) match(common_sids, sv))

    n_samp_out <- length(common_sids)

    bdose_con <- file(bdose_file, open = "wb")
    on.exit(close(bdose_con), add = TRUE)
    writeBin(.bdose5_magic, bdose_con)

    indices     <- numeric(sum(n_snps))
    current_pos <- 4.0
    snp_count   <- 0L

    for (k in seq_len(n_files)) {
      si <- samp_idx[[k]]
      for (i in seq_len(n_snps[k])) {
        snp <- getsnp(bdinfos[[k]], i, dosageonly = FALSE)
        ds  <- snp$dosage[si]
        if (all(is.na(snp$p0))) {
          gp <- rep(NA_real_, 3L * n_samp_out)
        } else {
          gp <- as.numeric(rbind(snp$p0[si], snp$p1[si], snp$p2[si]))
        }
        snp_count            <- snp_count + 1L
        block                <- compress_snp_block(ds, gp)
        indices[snp_count]   <- current_pos
        writeBin(block, bdose_con)
        current_pos <- current_pos + length(block)
      }
    }

    # Output SNPs: all SNPs from all files in file order
    out_snps <- do.call(rbind, lapply(bdinfos, function(b) b$snps))
    rownames(out_snps) <- NULL

    # Output samples: common subjects in first-file order
    out_samples <- bdinfos[[1L]]$samples[samp_idx[[1L]], , drop = FALSE]
    rownames(out_samples) <- NULL
  }

  # --- Save output bdinfo ---
  ref      <- bdinfos[[1L]]
  new_bdinfo <- list(
    filename       = normalizePath(bdose_file),
    usesfid        = ref$usesfid,
    samples        = out_samples,
    onechr         = length(unique(out_snps$chromosome)) == 1L,
    snpidformat    = ref$snpidformat,
    snps           = out_snps,
    snpinfo        = list(),
    additionalinfo = structure(
      list(format     = 5L,
           subformat  = 1L,
           headersize = 4L,
           numgroups  = 1L,
           groups     = n_samp_out),
      class = "bdose-info"
    ),
    datasize       = integer(0),
    indices        = indices
  )
  class(new_bdinfo) <- "genetic-info"
  saveRDS(new_bdinfo, paste0(bdose_file, ".bdi"))

  invisible(NULL)
}
