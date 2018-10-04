#' Randomize Headrow Trays
#' 
#' @param ents Vector listing names of entries to include in the experiment
#' @param nreps Integer specifying the number of replications 
#'   (i.e. full blocks)
#' @param tray_cells Integer specifying the number of cells contained by each
#'   headrow tray
#' @param tray_rows Integer specifying the number of rows per range. Note that
#'   this "rows" is shorthand for \emph{headrows}, so this might more 
#'   intuitively be thought of as columns.
#' @param start_cell Integer specifying the starting number for the first cell
#' @param randomize_first Logical indicating whether or not to randomize the
#'   first replication
#' @param checks A list with two elements:
#'   \enumerate{
#'     \item A string or character vector listing the names of check(s) to
#'     include
#'     \item An integer specifying how often to repeat the checks  
#'   }
#' @param seed Integer to seed the random number generator. Default NA value
#'   will not seed the random number generator.
#' @return A dataframe containing the randomized headrows and all pertinent
#'   identifying information
#' @details The headrows will be returned with randomized replications of
#'   experimental entries occuring one after another, with check entries inserted
#'   at regular intervals. Note that this returns a linear arrangement of trays
#'   with global range values increasing to infinity (i.e., the headrow trays
#'   are not wrapped at the length of the experiment field). Checks are first 
#'   inserted at the beginning (i.e. tray 1, range 1, row 1), and then repeated 
#'   at the specified number of cells. As noted above, the somewhat unfortunate
#'   nomenclature of head\emph{rows} means that what would intuitively be
#'   thought of as "rows" are labeled "ranges", and what would intuitively be
#'   thought of as "columns" are labeled "rows".
#' @export
randomize_headrows <- function(ents, nreps = 1, tray_cells = 80, tray_rows = 4,
                               start_cell = 1, randomize_first = TRUE,
                               checks = list("CHECK", 40), seed = NA) {

  if (!is.na(seed)) { set.seed(seed) }
  
  ## Sanity checks
  if (tray_cells %% tray_rows != 0) {
    stop("tray_cells must be a multiple of tray_rows")
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
  if (!randomize_first) { reps[[1]]$entry <- ents }
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
  if(nrow(headrows) %% tray_rows != 0) {
    nfill <- tray_rows - (nrow(headrows) %% tray_rows)
    fill_df <- data.frame("rep" = max(headrows$rep),
                          "entry_type" = "fill",
                          "entry" = rep("FILL", nfill))
    headrows <- rbind(headrows, fill_df)
  }
  
  ## Now add in additional identifying columns
  ntrays <- ceiling(nrow(headrows) / tray_cells)
  headrows$tray <- rep(1:ntrays, each = tray_cells)[1:nrow(headrows)]
  headrows$tray_cell <- rep_len(1:tray_cells, length.out = nrow(headrows))
  headrows$row <- rep_len(1:tray_rows, length.out = nrow(headrows))
  nrange <- tray_cells / tray_rows
  range_vec <- rep(1:nrange, each = tray_rows)
  headrows$tray_range <- rep_len(range_vec, length.out = nrow(headrows))
  headrows$exp_range <- rep(1:(ntrays*nrange), each = tray_rows)[1:nrow(headrows)]
  headrows$exp_cell <- seq(1:nrow(headrows)) + start_cell - 1
  
  ## Rearrange columns and return
  headrows <- headrows[c("rep", "tray", "tray_range", "row", "tray_cell", 
                         "exp_range", "exp_cell", "entry_type", "entry")]
  return(headrows)
}