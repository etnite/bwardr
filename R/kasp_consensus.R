#' Generate consensus calls for KASP marker data
#' 
#' @param kasp A dataframe containing KASP-formatted marker data. This should
#'   consist of at least one column with genotype names, and additional columns
#'   containing "X:X" and "Y:Y" for the two homozygous states, "X:Y" or "Y:X"
#'   for the heterozgous states, and "NO CALL" or "NOCALL" for missing data.
#' @param geno_col String identifying which column in kasp contains genotype names
#' @param data_cols Character vector identifying which columns in kasp will be
#'   retained and treated as marker data columns
#' @param ambiguous String consisting of either "strict" or "lax". If "strict",
#'   marker calls with uncertainty (followed by a question mark) will be converted
#'   to missing data. If "lax", uncertain calls will be converted to their
#'   putative call (e.g. "X:X?" will be converted to "X:X")
#' @param out_format String consisting of either "KASP", "VCF", or "numeric". If
#'   "KASP", output will be in same format as the input. IF "VCF", output will
#'   be encoded as "0/0" and "1/1" for the homozgous states, "0/1" for the
#'   heterozygous state, and "./." for missing data. If "numeric", output will be
#'   encoded as -1 and 1 for homozygous states, 0 for the heterozygous state,
#'   and NA for missing data.
#' @return Data frame with consensus calling performed, encoded as specified by
#'   out_format
#' @details kasp_consensus generates consensus calls for KASP marker data when
#'   some lines have been genotyped multiple times. It does this using a simple
#'   algorithm in which markers are first converted to numeric format, and then 
#'   subsequently summed for each genotype. Note that it cannot differentiate
#'   between randomly missing data and true null alleles. The input data must be
#'   reasonably well-formatted. The function can handle a variety of different
#'   formats (e.g. "x:x", "X:X", "NO CALL", "nocall"), as well as uncertain calls
#'   (e.g. "X:X?") but other more human-friendly data input will simply be
#'   changed to missing data, with a warning that data has been converted to NA
#'   by coercion.
#' @export
kasp_consensus <- function(kasp, geno_col, data_cols, 
                           ambiguous = "strict", out_format = "KASP") {
  
  ## Sanity checks
  stopifnot(out_format %in% c("KASP", "VCF", "numeric"))
  stopifnot(ambiguous %in% c("strict", "lax"))
  
  ## Split input data frame into genotypes column and data columns;
  ## convert data columns into vector
  genos <- kasp[, geno_col]
  k_dat <- unlist(kasp[, data_cols])
  k_dat <- toupper(k_dat)
  
  ## Either preserve ambiguous calls, or set to missing data
  ## depending on value of 'ambiguous'
  if (ambiguous == "lax") {
    k_dat <- gsub("?", "", k_dat, fixed = TRUE)
  } else if (ambiguous == "strict") {
    k_dat[grep("?", k_dat, fixed = TRUE)] <- NA
  }
  
  ## Convert to numeric format
  k_dat[k_dat == "X:X"] <- -1
  k_dat[k_dat == "Y:Y"] <- 1
  k_dat[k_dat == "X:Y" | k_dat == "Y:X"] <- 0
  k_dat[k_dat == "NO CALL" | k_dat == "NOCALL"] <- NA
  k_dat[k_dat == "NULL"] <- NA
  
  ## Convert back to data frame
  k_dat <- as.data.frame(matrix(k_dat, ncol = length(data_cols), dimnames = list(NULL, data_cols)))
  
  ## Combine line names and data columns - set data cols to integer class
  k_num <- data.frame("Line" = genos, k_dat, stringsAsFactors = FALSE)
  k_num[, data_cols] <- lapply(k_num[, data_cols], as.integer)
  
  ## Generate similar dataframe to record true hets
  hets <- k_num
  hets[, data_cols] <- abs(hets[, data_cols])
  
  ## In the main df, missing data and hets are set to 0
  ## In the hets df, only hets are set to 0
  k_num[is.na(k_num)] <- 0
  hets[is.na(hets)] <- 1
  
  ## Perform summation for each genotype
  k_sum <- stats::aggregate(. ~ Line, data = k_num, FUN = sum)
  hets_sum <- stats::aggregate(. ~ Line, data = hets, FUN = sum)
  
  ## Once again we split the df by genotype column and data cols
  genos <- k_sum$Line
  k_dat <- k_sum[, data_cols]
  
  ## Differentiate missing data and true hets.
  k_dat[k_dat == 0] <- NA
  for (i in data_cols) {
    k_dat[[i]][hets_sum[[i]] == 0] <- 0
  }
  k_dat[k_dat < 0] <- -1
  k_dat[k_dat > 0] <- 1
  
  ## Convert to selected output format
  ## If "numeric" is selected for out_format, then k_dat is unchanged
  if (out_format == "KASP") {
    k_dat[k_dat == -1] <- "X:X"
    k_dat[k_dat == 1] <- "Y:Y"
    k_dat[k_dat == 0] <- "X:Y"
    k_dat[is.na(k_dat)] <- "NO CALL"
  } else if (out_format == "VCF") {
    k_dat[k_dat == -1] <- "0/0"
    k_dat[k_dat == 1] <- "1/1"
    k_dat[k_dat == 0] <- "0/1"
    k_dat[is.na(k_dat)] <- "./."
  }
  
  ## Recombine to produce the final output df
  out_df <- data.frame(genos, k_dat, stringsAsFactors = FALSE)
  names(out_df)[names(out_df) == "genos"] <- geno_col
  return(out_df)
}