#' Convert Human-Readable DNA Length to Number of Nucleotides
#'
#' This function takes a human-readable DNA length and converts it to the actual
#' number of nucleotides. Supported units are 'Mb', 'kb', and 'bp'.
#' If no unit is provided, it assumes the unit is 'bp'.
#'
#' @param dna_length The human-readable DNA length (e.g., "0.5 Mb") or a numeric value.
#'
#' @return The DNA length in number of nucleotides.
#' @export
#'
#' @examples
#' \dontrun{
#' convert_dna_length("0.5 Mb") # Returns 500000
#' convert_dna_length("1000") # Returns 1000
#' convert_dna_length(1000) # Returns 1000
#' convert_dna_length("2 kb") # Returns 2000
#' }
convert_dna_length <- function(dna_length) {
  # If the input is numeric, assume the unit is 'bp'
  if (base::is.numeric(dna_length)) {
    return(dna_length)
  }

  # Remove leading and trailing whitespace
  dna_length <- base::gsub("^\\s+|\\s+$", "", dna_length)

  # Extract the numeric part and the unit
  if (base::grepl(" ", dna_length)) {
    parts <- base::strsplit(dna_length, " ")[[1]]
    num <- base::as.numeric(parts[1])
    unit <- parts[2]
  } else {
    # Handle case where no space exists between the number and the unit
    num <- base::gsub("[^0-9\\.]", "", dna_length) %>% base::as.numeric()
    unit <- base::gsub("[0-9\\.]", "", dna_length)
  }

  # If no unit is provided, assume 'bp'
  if (unit == "") {
    unit <- "bp"
  }

  # Check that the unit is one of the allowed units ('Mb', 'kb', 'bp')
  if (!unit %in% c("Mb", "kb", "bp")) {
    base::stop("Invalid unit. Allowed units are 'Mb', 'kb', and 'bp'.")
  }

  # Convert to the number of nucleotides
  if (unit == "Mb") {
    return(num * 1e6)
  } else if (unit == "kb") {
    return(num * 1e3)
  } else {
    return(num)
  }
}
