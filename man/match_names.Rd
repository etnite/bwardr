% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.R
\name{match_names}
\alias{match_names}
\title{Approximate matching between character vectors}
\usage{
match_names(a, b, max_dist = 5, drop_incomparables = TRUE, drop_nonsubs = TRUE)
}
\arguments{
\item{a}{A character vector}

\item{b}{A character vector}

\item{max_dist}{Integer - the maximum edit distance to consider for approximate
string matching using the longest common substring method}

\item{drop_incomparables}{Logical - whether to output names of entries in either
input character vector that could not be approximately matched to any entry
in the other input vector}

\item{drop_nonsubs}{Logical - whether to drop entries in either input character
vector that are not a full substring of any entry in the other input vector}
}
\value{
A dataframe containing seven columns: The original and simplified
  versions of the input vectors (see details), the edit distance between the
  two simplified vector entries, the difference between the two simplified
  vector entries, and a logical column indicating whether an entry in one
  of the vectors is a full substring of the corresponding entry in the other
}
\description{
Approximate matching between character vectors
}
\details{
For two input character vectors, the function will first create
  "simplified" versions of each, in which all non-alphanumeric characters are
  stripped out, and all alphabetic characters are converted to uppercase.
  Approximate string matching is then performed using the longest common
  substring (LCS) method in the stringdist package. Note that this may produce
  one-to-multiple matchings. The "substring" column in the output dataframe
  will indicate whether one term is a substring fully contained within the
  other. For instance, "Foobar" and "Foobaz" share the substring "Fooba",
  but neither is a full substring of the other. In contrast, "Foobar" is a valid
  substring of "Foobarbaz". The difference between any two corresponding
  elements in the input vectors is calculated using base R's adist() function.
  This uses the Levenshtein distance rather than the lcs metric, though its use
  here is "protected" by the previous use of the LCS metric to perform the
  approximate string matching.
}
