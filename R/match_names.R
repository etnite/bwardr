#' Approximate matching between character vectors
#' 
#' @param a A character vector
#' @param b A character vector
#' @param max_dist Integer - the maximum edit distance to consider for approximate
#'   string matching using the longest common substring method
#' @param drop_incomparables Logical - whether to output names of entries in either
#'   input character vector that could not be approximately matched to any entry
#'   in the other input vector
#' @param drop_nonsubs Logical - whether to drop entries in either input character
#'   vector that are not a full substring of any entry in the other input vector
#' @return A dataframe containing six columns: The original and simplified
#'   versions of the input vectors (see details), the edit distance between the
#'   two simplified vector entries, and a logical column indicating whether an entry in one
#'   of the vectors is a full substring of the corresponding entry in the other
#' @details For two character input vectors, the function will first create
#'   "simplified" versions of each, in which all non-alphanumeric characters are
#'   stripped out, and all alphabetic characters are converted to uppercase.
#'   Approximate string matching is then performed using the longest common
#'   substring (lcs) method in the stringdist package. Note that this may produce
#'   one-to-multiple matchings. The "substring" column in the output dataframe
#'   will indicate whether one term is a substring fully contained within the
#'   other. For instance, "Foobar" and "Foobaz" share the substring "Fooba",
#'   but neither is a full substring of the other. In contrast, "Foobar" is a valid
#'   substring of "Foobarbaz".
#' @export
match_names <- function(a, b, max_dist = 5, drop_incomparables = TRUE, drop_nonsubs = TRUE) {
  
  ## Perform name simplification for the a vector
  a_df <- data.frame("a_original" = a)
  a_df$a_simple <- gsub("[^A-Za-z0-9]", "", toupper(a_df$a_original))
  
  ## Perform name simplification for the b vector
  b_df <- data.frame("b_original" = b)
  b_df$b_simple <- gsub("[^A-Za-z0-9]", "", toupper(b_df$b_original))
  
  ## Get approximate matches and edit distances
  a_df$b_simple <- b_df$b_simple[stringdist::amatch(a_df$a_simple, b_df$b_simple, method = "lcs", maxDist = max_dist)]
  a_df$ab_dist <- stringdist::stringdist(a_df$a_simple, a_df$b_simple, method = "lcs")
  
  ## Find which matched strings are valid, full substrings of one another
  a_in_b <- b_in_a <- rep(NA, nrow(a_df))
  for (i in 1:nrow(a_df)) {
    a_in_b[i] <- grepl(a_df$a_simple[i], a_df$b_simple[i])
    b_in_a[i] <- grepl(a_df$b_simple[i], a_df$a_simple[i])
  }
  a_df$substring <- as.logical(a_in_b + b_in_a)
  
  if (drop_incomparables) {
    a_df <- merge(a_df, b_df, by = "b_simple", all = FALSE)
  } else {
    a_df <- merge(a_df, b_df, by = "b_simple", all = TRUE)
  }
  
  ## Ensure substring col. only contains TRUE/FALSE, no NA values
  a_df$substring[is.na(a_df$substring)] <- FALSE
  
  if (drop_nonsubs) {
    a_df <- a_df[a_df$substring, ]
  }
  
  return(a_df[c("a_original", "a_simple", "b_original", "b_simple", "ab_dist", "substring")])
}