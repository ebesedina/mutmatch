#' Shuffle Mutation Counts for a Single Cohort
#'
#' This function takes a mutation table for a single cohort and returns a table with permuted mutation counts.
#' The shuffling is done based on the marginal sums of mutations and ntAtRisk for each mutation type.
#' The returned table will have a new column named 'permutation_id' to identify each permutation.
#'
#' @param mutation_table A data.table containing mutation information. It should contain only one category (cohort)
#' except for the 'isTarget' and 'Mutation' variables.
#' @param ntimes The number of times the shuffling is to be done.
#'
#' @return A data.table containing the permuted mutation counts.
#' @export
#'
#' @examples
#' \dontrun{
#' shuffled_data <- shuffle_mutation_counts(mutation_table, ntimes = 100)
#' }
shuffle_mutation_counts <- function(mutation_table, ntimes = 50) {
  # Set seed for reproducibility
  set.seed(42)

  # Calculate marginal sums of mutations and ntAtRisk for each mutation type
  mutation_table_marginal_sums <- mutation_table %>%
    dplyr::group_by(Mutation) %>%
    dplyr::summarize(
      MutationNumber = sum(MutationNumber),
      ntAtRisk = sum(ntAtRisk)
    ) %>%
    data.table::as.data.table()

  # Identify unique mutation types present in the data
  mutation_types_present <- unique(mutation_table_marginal_sums$Mutation)

  # Set keys for fast subsetting
  data.table::setkey(mutation_table, Mutation)
  data.table::setkey(mutation_table_marginal_sums, Mutation)

  # Perform the shuffling
  mutation_table_permutated <- data.table::rbindlist(lapply(mutation_types_present, function(mutation_type) {
    original_table <- mutation_table[.(mutation_type)]
    original_table <- original_table[order(isTarget)]
    data.table::setkey(original_table)

    if (!length(unique(original_table$isTarget)) == 2) {
      permutated_table_n <- original_table[rep_len(
        seq_len(nrow(original_table)),
        nrow(original_table) * ntimes
      )]
      permutated_table_n[, permutation_id := rep(seq_len(ntimes), each = nrow(original_table))]
    } else {
      marginal_sums <- mutation_table_marginal_sums[.(mutation_type)]
      marginal_sum_MutationNumber <- marginal_sums$MutationNumber
      marginal_sum_ntAtRisk <- marginal_sums$ntAtRisk

      central_gene_counts <- stats::rbinom(
        n = ntimes,
        size = marginal_sum_MutationNumber,
        prob = original_table[isTarget == 1]$ntAtRisk / marginal_sum_ntAtRisk
      )

      # Get the number of rows in the original table, should be 2
      num_rows_original <- nrow(original_table)

      # Replicate the entire original_table ntimes
      permutated_table_n <- original_table[rep(seq_len(.N), ntimes)]

      # Add a permutation ID
      permutated_table_n[, permutation_id := rep(1:ntimes, each = num_rows_original)]

      # Update MutationNumber for isTarget == 1 and isTarget == 0
      permutated_table_n[isTarget == 1, MutationNumber := central_gene_counts]
      permutated_table_n[isTarget == 0, MutationNumber := marginal_sum_MutationNumber - central_gene_counts]
    }

    data.table::setcolorder(permutated_table_n, c(colnames(mutation_table), "permutation_id"))
    return(permutated_table_n)
  }))

  # Convert 'isTarget' and 'Mutation' to factors
  mutation_table_permutated[, isTarget := as.factor(isTarget)]
  mutation_table_permutated[, Mutation := as.factor(Mutation)]
  mutation_table_permutated[]
  return(mutation_table_permutated)
}
