#' Clean a vector of delimited strings
#'
#' @param vec A vector of delimited strings to clean
#' @param sep String - delimiter used within each element of the vector to be cleaned
#' @return Vector - the cleaned input vector, with each element consisting of 
#'   sorted, unique elements
#' @details This function handles cases where delimited data is stored within
#'   each element of a vector. Such vectors are in a sense vectors of
#'   vectors, though each element is stored as a delimited string. In addition, these
#'   vectors resemble lists, which can store a sequence of vectors. The function
#'   will loop through each element of the input vector, check for the presence of a
#'   string containing the specified delimiter, split the element string apart by
#'   the delimiter, attempt to convert to numeric, and then sort the unique values
#'   contained within the split string.
#' @examples
#' df <- data.frame("mycol" = c("B; A; M", "D; O; G", "6; 2; 9; 10; 12"), 
#'                  stringsAsFactors = FALSE)
#' head(df$mycol)
#' df$mycol <- clean_dvec(df$mycol, sep = "; ")  
#' head(df$mycol)               
#' @export
clean_dvec <- function(vec, sep = ";") {
  stopifnot(is.vector(vec))
  for (i in 1:length(vec)) {
    if (!is.na(vec[i]) & grepl(sep, vec[i])) {
      element_vec <- strsplit(vec[i], split = sep, fixed = TRUE)[[1]]
      
      ## Try to convert to numeric; if NAs are introduced, then abort
      num_conv <- suppressWarnings(as.numeric(element_vec))
      if (sum(is.na(num_conv)) == sum(is.na(element_vec))) {
        element_vec <- num_conv
      }
      
      ## Collapse into string after sorting unique values
      vec[i] <- paste(sort(unique(element_vec)), collapse = sep)
    }
  }
  return(vec)
}
