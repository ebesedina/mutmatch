#' Generate Possible Mutation Spectra
#'
#' This function generates all possible mutation spectra given a nucleotide context length.
#' It constructs trinucleotide or pentanucleotide combinations, expands them with possible mutations,
#' and then filters the results to create a list of possible mutations.
#'
#' @param nt_length A character string specifying the nucleotide context length for constructing the mutation spectra.
#' Options are "tri" for trinucleotides (MS96) and "penta" for pentanucleotides. Default is c("tri", "penta").
#'
#' @return A data table with columns "Context" and "Mutation". The "Context" column contains the nucleotide
#' context in which the mutation occurs, and the "Mutation" column contains the mutation represented as "Context>MutatedNucleotide".
#'
#' @examples
#' \dontrun{
#' get_mutation_spectra() # By default, uses trinucleotide context
#' get_mutation_spectra(nt_length = "penta") # Uses pentanucleotide context
#' }
#' @export
get_mutation_spectra <- function(nt_length = c("tri", "penta")) {
  nt_length <- base::match.arg(nt_length)

  # Generate combinations based on the oligonucleotide length
  if (nt_length == "tri") {
    nt_length <- 3
    pos_comparison <- c(2, 5)
  } else if (nt_length == "penta") {
    nt_length <- 5
    pos_comparison <- c(3, 7)
  } else {
    stop("The 'nt_length' parameter should be either 'tri', or 'penta'")
  }

  # Define the nucleotide alphabet and initialize a data frame for possible mutations
  nucleotide_alphabet <- c("A", "C", "G", "T")

  # Generate pentanucleotide combinations using CJ
  oligonucleotide_combinations <- do.call(data.table::CJ, replicate(nt_length, nucleotide_alphabet, simplify = FALSE))
  oligonucleotide_strings <- do.call(paste, c(oligonucleotide_combinations, sep = ""))

  # Expand to include the nucleotide alphabet and collapse to form mutations
  possible_mutations <- data.table::CJ(oligonucleotide = oligonucleotide_strings, mutation = nucleotide_alphabet)
  possible_mutations[, mutation_spectra := paste(oligonucleotide, mutation, sep = ">")]

  # Filter based on the original nucleotide and the mutatated one
  original_nucleotide <- substr(possible_mutations$mutation_spectra, pos_comparison[1], pos_comparison[1])
  mutated_nucleotide <- substr(possible_mutations$mutation_spectra, pos_comparison[2], pos_comparison[2])
  condition <- (original_nucleotide != mutated_nucleotide) & !(original_nucleotide %in% c("A", "G"))

  # Create the final data frame
  possible_mutations <- data.table::data.table(
    "Context" = possible_mutations$oligonucleotide[condition],
    "Mutation" = possible_mutations$mutation_spectra[condition]
  )
  return(possible_mutations)
}
