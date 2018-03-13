#' Generate GxE Plot for a trait
#'
#' @param trait_df - a dataframe or tibble with at least three columns - 
#'   GENOTYPE, ENV, and one or more columns for trait values
#' @param trait_of_int - string matching the column name for the trait of
#'   interest in trait_df
#' @param min_core - integer specifying the minimum number of "core" genotypes
#'   (i.e. genotypes tested across all environments) required to produce plot. If
#'   the number of core genotypes falls below this number, NULL will be returned
#' @return a ggplot plot object, or NULL
#' @details This function generates GxE plots in a manner similar to that of
#'   Finlay and Wilkinson (1963; https://doi.org/10.1071/AR9630742). Briefly,
#'   a "core" set of genotypes tested across all environments is identified.
#'   The environments are then arranged in increasing order by the mean values
#'   of the core genotypes within each environment, and then a regression is fit
#'   for each individual core genotype.
#' @importFrom magrittr %>%
#' @export
ge_plot_gen <- function(trait_df, trait_of_int, min_core = 5) {
  
  ## Identify and subset core lines
  genos_by_env <- list()
  for (strat in unique(trait_df$ENV)) {
    genos_by_env[[strat]] <- trait_df$GENOTYPE[trait_df$ENV == strat]
  }
  core_lines <- Reduce(intersect, genos_by_env)
  
  if (length(core_lines) >= min_core) {
    dplyr::filter(trait_df, GENOTYPE %in% core_lines) %>%
      dplyr::select(one_of(c("GENOTYPE", "ENV", trait_of_int))) ->
      sub_df
    
    ## Here standard evaluation is used (!!sym())
    dplyr::group_by(sub_df, ENV) %>%
      dplyr::summarise(Mean = mean(!!sym(trait_of_int))) %>%
      dplyr::arrange(Mean) %>%
      dplyr::mutate(GENOTYPE = "MEAN") ->
      env_order
    
    ## Add the within-env means to the data
    colnames(env_order)[colnames(env_order) == "Mean"] <- trait_of_int
    sub_df <- dplyr::bind_rows(env_order, sub_df)
    sub_df$type <- "entry"
    sub_df$type[sub_df$GENOTYPE == "MEAN"] <- "mean"
    
    ## Order the environments by their means
    env_order <- as.character(env_order$ENV)
    sub_df$ENV <- factor(sub_df$ENV, levels = env_order)
    sub_df <- sub_df[order(sub_df$ENV), ]
    
    ## Generate GxE plot
    ggplot2::ggplot(sub_df, aes_string(x = "ENV", y = trait_of_int)) +
      geom_smooth(aes(group = GENOTYPE, color = type), method = "lm", size = 0.75, se = FALSE) +
      scale_color_manual(values = c(mean = "red", entry = "black")) +
      xlab("Environment") +
      ylab(trait_of_int) +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, vjust = 0.25)) ->
      ge_plot
    
    return(ge_plot)
    
  } else {
      warning(paste("Number of core genotypes falls below specified threshold for trait", trait_of_int))
      return(NULL)
  }
}