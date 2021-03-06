#' Standardize character vectors to a common format
#' 
#' @param str_in Character vector
#' @param sep Character defining the word delimiter for the output character
#'   vector, either " " (space), ".", "-", "_", or "" (no delimiter)
#' @param case String consisting of either "upper" or "lower" to select whether
#'   to convert output to all uppercase or all lowercase
#' @return character vector of same dimensions as input vector, but with all
#'   elements updated to the new format
#' @details This function somewhat crudely attempts to force character vectors
#'   to a common format. It does this by converting everything to upper or lowercase, 
#'   and replacing all common word delimiters (space, period, underscore, dash, 
#'   comma) with the selected one.
#'   Note that of the three, only underscores or periods will form syntactically-
#'   valid R variable names (override by reading in tables with check.names = FALSE). 
#'   However, underscores can cause problems with reading VCFs in PLINK, as they 
#'   will by default be interpreted as delimiters for family and individual IDs.
#'   Choosing spaces or commas as the output delimiter is generally not recommended, as this
#'   can cause issues in whitespace-delimited or comma-delimited text files, respectively.
#' @export
stand_str <- function(str_in, sep = "-", case = "upper") {
  
  ## Check for appropriate word delimiter selection
  if(!sep %in% c(" ", "-", "_", ".", "", ",")) {
    stop('Please select one of the following for the output word delimiter:
         " " [space], ".", "-", "_", ",", "" [no delimiter]')
  }
  
  ## Remove all non-ASCII characters
  str_out <- iconv(str_in, "latin1", "ASCII", sub="")
  
  ## Convert case
  if (case == "upper") {
      str_out <- toupper(str_out)
  } else if (case == "lower") {
      str_out <- tolower(str_out)
  }
  
  ## Remove leading/trailing whitespace
  str_out <- trimws(str_out, which = "both")
  
  ## Remove comment characters (hash/pound sign)
  str_out <- gsub("#", "", str_out, fixed = TRUE)
  
  ## Convert newlines to dashes
  str_out <- gsub("[\r\n]", "-", str_out)
  
  ## Convert all delimiters (except no delimiter) to dashes
  str_out <- gsub("\\s|-|,|\\.|_", "-", str_out)
  
  ## Replace repeated dashes with single
  str_out <- gsub("-+", "-", str_out)
  
  ## Remove leading/trailing dashes
  str_out <- gsub("^-|-$", "", str_out)
  
  ## Convert dashes to final delimiter choice
  str_out <- gsub("\\s|-|\\.|_", sep, str_out)
  
  return(str_out)
}
