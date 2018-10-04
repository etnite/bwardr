#' Randomize Headrow Trays
#' 
#' @param ents A vector listing names of entries to include in the experiment
#' @param nreps An integer specifying the number of replications 
#'   (i.e. full blocks)
#' @param tray_cells An integer specifying the number of cells contained by each
#'   headrow tray
#' @param tray_cols An integer specifying the number of columns in a headrow
#'   tray
#' @param checks A list with two elements:
#'   \enumerate{
#'     \item A string or character vector listing the names of check(s) to
#'     include
#'     \item An integer specifying how often to repeat the checks  
#'   }
#' @param seed An optional integer to seed the random number generator
#' @return A dataframe containing the randomized headrows and all pertinent
#'   identifying information
#' @details The headrows will be returned with randomized replications of
#'   experimental entries occuring one after another, with check entries inserted
#'   at regular intervals. Note that this returns a linear arrangement of trays
#'   with global range values increasing to infinity (i.e., the headrow trays
#'   are not wrapped at the length of the experiment field). Checks are first 
#'   inserted at the beginning (i.e. tray 1, range 1, col 1), and then repeated 
#'   at the specified number of cells.
#' @export
randomize_headrows <- function(ents, nreps, tray_cells = 80, tray_cols = 4,
                               checks = list("CHECK", 40), seed = NA) {

  if (!is.na(seed)) { set.seed(seed) }
  
  ## Sanity checks
  if (tray_cells %% tray_cols != 0) {
    stop("tray_cells must be a multiple of tray_cols")
  }
  
  input_ents <- ents
  ents <- sort(unique(ents))
  if (length(ents) != length(input_ents)) {
    warning("Duplicated entries were removed from 'ents' input vector")
  }
  
  ## Create dataframe of concatenated randomized replications
  reps <- vector("list", length = nreps)
  for(i in 1:nreps) {
    reps[[i]] <- data.frame("rep" = i,
                            "entry_type" = "experimental",
                            "entry" = sample(ents, size = length(ents), replace = FALSE))
  }
  reps_df <- do.call("rbind", reps)
  
  ## This snippet from StackExchange 
  ## https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows?lq=1
  ## split a dataframe into equal sized chunks (with "remainder")
  chunk <- checks[[2]] - length(checks[[1]])
  n <- nrow(reps_df)
  r  <- rep(1:ceiling(n/chunk), each=chunk)[1:n]
  exp_list <- split(reps_df, r)
  
  ## Create the dataframe of checks to insert into the entries dataframe
  checks_df <- data.frame("rep" = NA,
                          "entry_type" = "check",
                          "entry" = checks[[1]])
  
  ## Loop through the list of entries chunks, and insert the checks between them
  for (i in 1:length(exp_list)) {
    ents_rep <- exp_list[[i]]$rep[1]
    check_insert <- checks_df
    check_insert$rep <- ents_rep
    exp_list[[i]] <- rbind(check_insert, exp_list[[i]])
  }
  headrows <- do.call("rbind", exp_list)
  
  ## Add in fill if a partial range remains at the end
  if(nrow(headrows) %% tray_cols != 0) {
    nfill <- tray_cols - (nrow(headrows) %% tray_cols)
    fill_df <- data.frame("rep" = max(headrows$rep),
                          "entry_type" = "fill",
                          "entry" = rep("FILL", nfill))
    headrows <- rbind(headrows, fill_df)
  }
  
  ## Now add in additional identifying columns
  ntrays <- ceiling(nrow(headrows) / tray_cells)
  headrows$tray <- rep(1:ntrays, each = tray_cells)[1:nrow(headrows)]
  headrows$cell <- rep_len(1:tray_cells, length.out = nrow(headrows))
  headrows$col <- rep_len(1:tray_cols, length.out = nrow(headrows))
  nrange <- tray_cells / tray_cols
  range_vec <- rep(1:nrange, each = tray_cols)
  headrows$range <- rep_len(range_vec, length.out = nrow(headrows))
  headrows$global_range <- rep(1:(ntrays*nrange), each = tray_cols)[1:nrow(headrows)]
  
  ## Rearrange columns and return
  headrows <- headrows[c("rep", "tray", "range", "global_range", "col", "cell", 
                         "entry_type", "entry")]
  return(headrows)
}