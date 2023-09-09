#' Retrieve mutations within a genomic range
#'
#' This function identifies mutations within a given genomic range (grObject)
#' and stratifies the result based on the annotation and cluster number.
#'
#' @inheritParams get_mutations_annot
#' @param annotation Data table with annotation details for each sample.
#'
#' @return A data table containing the mutations within the genomic range,
#' stratified based on annotation and cluster number.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume grObject, grMutations, and annotation are properly defined
#' result <- get_mutations_table(grObject, grMutations, annotation)
#' }
get_mutations_table <- function(grObject, grMutations, annotation = NULL) {
  # Mutations are stratified using MS96 (trinucleotide context)
  nt_length <- "tri"

  # Count the number of samples by cancer type
  sample_cohorts_number <- plyr::count(annotation[, -"Sample"])
  colnames(sample_cohorts_number)[colnames(sample_cohorts_number) == "freq"] <- "Sample_number"

  # Define mutation spectra
  possible_mutations <- get_mutation_spectra(nt_length = nt_length)
  annotation_cohorts <- annotation[, -"Sample"] %>% unique()
  possible_mutations_by_cohort <- dplyr::cross_join(annotation_cohorts, possible_mutations)

  # Calculate nucleotide composition
  oligonucleotide_composition <- calculate_trinucleotide_composition(grObject)

  # Define the cohorts and stratification features
  sample_cohorts <- base::setdiff(colnames(annotation), "Sample")
  strat_features <- base::c("Mutation", sample_cohorts)

  # Collect and Stratify Mutations
  mutationsData <- get_mutations_annot(grObject = grObject, grMutations = grMutations, nt_length = nt_length)
  mutationsData <- merge(mutationsData, annotation, by = "Sample")

  # Calculate mutation number for each cohort
  mutation_number_by_cohort <- plyr::count(mutationsData[, ..strat_features]) %>% data.table::data.table()
  colnames(mutation_number_by_cohort)[colnames(mutation_number_by_cohort) == "freq"] <- "MutationNumber"
  # selected_columns = c(strat_features, "Context")
  mutationsData <- base::merge(mutationsData[, ..strat_features] %>% unique(),
    mutation_number_by_cohort,
    by = strat_features
  )

  mutationsData <- base::merge(mutationsData,
    possible_mutations_by_cohort,
    by = strat_features, all.y = T
  )
  mutationsData[is.na(MutationNumber), MutationNumber := 0]

  # Calculate maximal mutation number for each cohort
  mutationsData <- base::merge(mutationsData, sample_cohorts_number, by = setdiff(strat_features, "Mutation"))
  mutationsData <- base::merge(mutationsData, oligonucleotide_composition[, c("Context", "Counts")], by = "Context")
  mutationsData[, ntAtRisk := Sample_number * Counts]

  # Convert stratification columns to factors and set "0" as the base level if it exists in the levels of the factor
  mutationsData[, (strat_features) := base::lapply(.SD, function(x) {
    x_factor <- as.factor(x)
    if ("0" %in% levels(x_factor)) {
      x_factor <- relevel(x_factor, ref = "0")
    }
    return(x_factor)
  }), .SDcols = strat_features]
  mutationsData <- mutationsData[, c(strat_features, "ntAtRisk", "MutationNumber"), with = F]

  return(mutationsData)
}
