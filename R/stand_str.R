#' Standardize character vectors to a common format
#' 
#' @param str_in Character vector
#' @param word_delim Character defining the word delimiter for the output character
#'   vector, either " " (space), ".", "-", or "_"
#' @return character vector of same dimensions as input vector, but with all
#'   elements updated to the new format
#' @details This function somewhat crudely attempts to force character vectors
#'   to a common format. It does this by converting everything to uppercase, and replacing
#'   all common word delimiters (space, period, underscore, dash) with the selected one.
#'   Note that of the three, only underscores or periods will form syntactically-
#'   valid R variable names (override by reading in tables with check.names = FALSE). 
#'   However, underscores can cause problems with reading
#'   VCFs in PLINK, as they will by default be interpreted as delimiters for
#'   family and individual IDs
#' @suggests stringdist, plyr
stand_str <- function(str_in, word_delim = "-") {
  
  ## Check for appropriate word delimiter selection
  if(!word_delim %in% c(" ", "-", "_", ".")) {
    stop('Please select one of the following for the output word delimiter:
         " " [space], ".", "-", "_"')
  }
  
  ## Convert everything to uppercase
  str_out <- toupper(str_in)
  
  ## Remove leading/trailing whitespace
  str_out <- trimws(str_out, which = "both")
  
  ## Replace emdashes with hyphens
  str_out <- gsub("--", "-", str_out, fixed = TRUE)
  
  ## Remove comment characters
  str_out <- gsub("#", "", str_out, fixed = TRUE)
  
  ## Convert asterices to dashes
  str_out <- gsub("*", "-", str_out, fixed = TRUE)
  
  ## Remove non-ASCII characters
  str_out <- iconv(str_out, "latin1", "ASCII", sub="")
  
  ## Swap out whitespace, dash, period, underscore for selected delimiter
  str_out <- gsub("\\s|-|\\.|_", word_delim, str_out)
  
  return(str_out)
}