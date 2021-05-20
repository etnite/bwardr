#' Convert between T3 input and output formats
#' 
#' @param df A dataframe containing either phenotypic data for import into the
#'   T3 database (see details), or phenotypic data exported from T3
#' @param direction String consisting of either "io" indicating that a dataframe
#'   in T3 input format has been supplied, and a dataframe in T3 export format
#'   should be supplied, or "oi" for the reverse conversion.
#' @details The T3 database \url{https://wheat.triticeaetoolbox.org} is a public
#'   Breedbase (\url{https://breedbase.org}) instance focusing on data collected
#'   from wheat, barley, and oats. T3 uses upload templates 
#'   (explained at \url{https://notes.triticeaetoolbox.org/s/h4qRqAfhe#trials-})
#'   for loading new data into the database. However, the data that is exported
#'   from the database has different column headers. This function is intended
#'   to convert between the input and output formats of the database. In lieu of
#'   any sort of dictionary object type in R, here we use a dataframe of key-value
#'   pairs to perform the substitution. This dataframe is hardcoded into the
#'   function. Note that not all columns in the T3 input format are represented
#'   in the output format. For more general key-value swapping in R, a function 
#'   such as plyr::mapvalues() is useful and much more performant than this function.
#' @return A list containing the following elements:
#' * df The updated dataframe
#' * notpresent A vector containing the keys (column names of the input dataframe)
#'   that have no corresponding values (column names to translate to in the output
#'   dataframe)
#' @export
#' @md
t3_in_out <- function(df, direction = "io") {
  
  stopifnot(direction %in% c("io", "oi"))
  
  ## Create a key-value dataframe for converting between T3 input/export formats
  keyval <- data.frame(
    "key" = c("trial_name", "breeding_program", "location", "year", 
               "design_type", "description", "trial_type", "plot_width", 
               "plot_length", "field_size","planting_date", "harvest_date", 
               "plot_name", "accession_name", "plot_number", "block_number", 
               "is_a_control", "rep_number", "range_number", "row_number", 
               "col_number", "seedlot_name", "num_seed_per_plot", "weight_gram_seed_per_plot",
               rep(NA, 17)),
    "val" = c("studyName", "programName", "locationName", "studyYear", "studyDesign", 
             "studyDescription", NA, "plotWidth", "plotLength", "fieldSize", 
             "plantingDate", "harvestDate", "observationUnitName", "germplasmName", 
             "plotNumber", "blockNumber", "entryType", "replicate", NA, 
             "rowNumber", "colNumber", "plantedSeedlotStockUniquename",
             "plantedSeedlotTransactionCount", "plantedSeedlotTransactionWeight", 
             "programDbId", "programDescription", "studyDbId", "fieldTrialIsPlannedToBeGenotyped", 
             "fieldTrialIsPlannedToCross", "locationDbId", "germplasmDbId", 
             "germplasmSynonyms", "observationLevel", "observationUnitDbId", 
             "plantNumber", "plantedSeedlotStockDbId", "plantedSeedlotCurrentCount", 
             "plantedSeedlotCurrentWeightGram", "plantedSeedlotBoxName", 
             "plantedSeedlotTransactionDescription", 
             "availableGermplasmSeedlotUniquenames")
  )

  ## Shape the keyval dataframe. If direction is "oi" then we switch the val/key cols
  ## Drop all keys that don't have a corresponding value
  if (direction == "oi") { names(keyval) <- c("val", "key") }
  keyval <- keyval[!is.na(keyval$val), ]
  keyval <- keyval[!is.na(keyval$key), ]
  
  ## Perform key/val replacement, or concatenate key to notpresent vector
  notpresent <- NULL
  for (i in 1:ncol(df)) {
    k <- names(df)[i]
    if (k %in% keyval$key) {
      names(df)[i] <- keyval$val[keyval$key == k]
    } else {
      notpresent <- c(notpresent, k)
    }
  }

  return(list("df" = df, "notpresent" = notpresent))
}