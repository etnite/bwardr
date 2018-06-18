#' Create common synonyms of strings
#'
#' @param s Character vector or string
#' @return Character vector containing common synonyms of the input string or
#'   character vector
#' @details This function takes an input string or vector of strings, 
#'   and creates an output vector that has uppercase, lowercase, and title case 
#'   versions of the input string(s) with different common delimiters 
#'   (space, hyphen, underscore, period). Note that the different strings present 
#'   in the output will have all delimiters the same (e.g. robin_sparrow_swift, 
#'   not robin.sparrow_swift)
#' @export
synonymize <- function(s) {
  
  ## Convert common delimiters to spaces, create vector of input string with
  ## lowercase, uppercase, and title case versions
  s <- tolower(s)
  s <- chartr("-_.", "   ", s)
  hi_lo <- c(s, toupper(s), stringr::str_to_title(s))
  
  ## Loop through the hi_lo vector, substituting in common delimiters
  delims <- c(" ", "-", "_", ".")
  delim_sub <- list()
  for (i in 1:length(delims)) {
    delim_sub[[i]] <- gsub(" ", delims[i], hi_lo, fixed = TRUE)
  }
  
  ## Concatenate everything into single vector and output
  out_vec <- unique(do.call("c", delim_sub))
  return(out_vec)
}
