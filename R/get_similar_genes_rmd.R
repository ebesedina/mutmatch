#' Load and Filter RmdGenesSimilarity Data
#'
#' This function loads chunks of RmdGenesSimilarity data from RDA files,
#' combines them, and then filters the combined data based on an
#' outlierNeighborsThreshold value.
#'
#' @param outlierNeighborsThreshold A numeric value used as the threshold for
#'        filtering the \code{logOR} column in the data. The absolute value of
#'        the \code{logOR} should be less than or equal to this threshold.
#'
#' @return A combined and filtered data.table object.
#'
#' @examples
#' # Example usage (assuming the function and RDA files exist)
#' # result_data <- .load_and_filter_RmdGenesSimilarity(0.25)
#'
.load_and_filter_RmdGenesSimilarity <- function(outlierNeighborsThreshold) {
  # Define the max thresholds for the chunks
  available_thresholds <- seq(0.1, 5.4, by = 0.1)

  # Determine which files need to be loaded
  thresholds_to_load <- available_thresholds[available_thresholds <= ceiling(outlierNeighborsThreshold * 10) / 10]

  # Reusable environment
  tmp_env <- new.env()

  # Use lapply to load data chunks
  loaded_data <- lapply(thresholds_to_load, function(threshold) {
    file_name <- paste0("RmdGenesSimilarity_max_", threshold, ".rda")
    full_path <- system.file("extdata", "RmdGenesSimilarity", file_name, package = "mutmatch")

    # Load the .rda file into the temporary environment
    load(full_path, envir = tmp_env)

    # Get the data chunk from the temporary environment
    chunk_name <- ls(envir = tmp_env)[1]
    tmp_env[[chunk_name]]
  })

  # Combine the loaded data chunks
  combined_data <- do.call(rbind, loaded_data)

  # Filter the combined data
  filtered_data <- combined_data[abs(combined_data$logOR) <= outlierNeighborsThreshold, ]

  return(filtered_data)
}


#' Get Similar Genes Based on Rmd (Regional Mutation Density)
#'
#' This function retrieves genes that are similar to a specified gene
#' (`hgnc`) based on a provided threshold (`outlierNeighborsThreshold`).
#' The function filters the data to return only the GeneIDs that are similar
#' to the specified gene based on the \code{logOR} value.
#'
#' @inheritParams get_gene_neighbors
#'
#' @return  character vector of unique GeneIDs that are similar to the specified gene.
#'
#' @examples
#' # Example usage (assuming the function and relevant data exist)
#' similar_genes <- get_similar_genes_rmd("KRAS", 0.2)
#' @export
get_similar_genes_rmd <- function(hgnc, outlierNeighborsThreshold) {
  # Load only the rows relevant to the specified outlierNeighborsThreshold to improve performance
  mutmatch_RmdGenesSimilarity <- .load_and_filter_RmdGenesSimilarity(outlierNeighborsThreshold)

  # Filter the data to find similar genes based on the outlierNeighborsThreshold and the specified gene (hgnc)
  genes_sim_rmd <- mutmatch_RmdGenesSimilarity %>%
    dplyr::filter((GeneName1 == hgnc | GeneName2 == hgnc) &
      abs(logOR) <= outlierNeighborsThreshold &
      GeneName1 != GeneName2) %>%
    dplyr::select(GeneID1, GeneID2) %>%
    base::unlist() %>%
    base::unique()

  return(genes_sim_rmd)
}
