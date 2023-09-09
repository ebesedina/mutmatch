#' Get Ensembl ID from HGNC Symbol
#'
#' This function retrieves the Ensembl ID corresponding to a given HGNC symbol.
#' It first attempts to find this in the `mutmatch_topTranscripts` data frame.
#' If unsuccessful, it checks for approved and alias symbols.
#'
#' @param hgnc Character string. The HGNC symbol for which the Ensembl ID is to be found.
#'
#' @return A string containing the Ensembl ID, or NULL if it cannot be found.
#'
#' @export
get_ensembl_from_hgnc <- function(hgnc) {
  # Try to find the Ensembl ID using the provided HGNC symbol
  ensembl_id <- mutmatch_topTranscripts %>%
    dplyr::filter(Hugo_Symbol == hgnc) %>%
    dplyr::pull(Gene) %>%
    base::unique()

  # If an Ensembl ID is found, return it
  if (length(ensembl_id) != 0) {
    return(ensembl_id)
  }

  # If no Ensembl ID is found, first check for approved symbol
  warning("Not found top-expressed transcripts for provided gene name, trying approved gene symbol instead.")
  gene_symbol_approved <- check_gene_aliases(hgnc, approved = TRUE)

  if (length(gene_symbol_approved) > 0) {
    ensembl_id <- mutmatch_topTranscripts %>%
      dplyr::filter(Hugo_Symbol %in% gene_symbol_approved) %>%
      dplyr::pull(Gene) %>%
      base::unique() %>%
      dplyr::first()
  }

  # If an Ensembl ID is found with an approved symbol, return it
  if (length(ensembl_id) == 1) {
    warning(paste("Continuing with", ensembl_id, "gene id for", hgnc))
    return(ensembl_id)
  } else {
    # If still not found, try searching for aliases
    warning("Not found top-expressed transcripts for approved symbols, trying possible gene aliases instead.")
    gene_symbol_aliases <- check_gene_aliases(hgnc, approved = FALSE)

    if (length(gene_symbol_aliases) > 0) {
      ensembl_id <- mutmatch_topTranscripts %>%
        dplyr::filter(Hugo_Symbol %in% gene_symbol_aliases) %>%
        dplyr::pull(Gene) %>%
        base::unique() %>%
        dplyr::first()
    }

    # If an Ensembl ID is finally found, return it
    if (length(ensembl_id) == 1) {
      warning(paste("Continuing with", ensembl_id, "gene id for", hgnc))
      return(ensembl_id)
    } else {
      # If Ensembl ID still can't be found, throw an error
      stop(paste("Not possible to find gene id with data about top-expressed transcripts for provided gene name", hgnc))
      return(NULL)
    }
  }
}

#' Check for Gene Aliases
#'
#' This function checks for aliases of a given gene symbol (HGNC) in the `mutmatch_geneAliases` data frame.
#' The function can return either approved aliases only or include previous and alias symbols.
#'
#' @param hgnc Character string. The HGNC symbol for which aliases are to be checked.
#' @param approved Logical. If TRUE, returns only the approved aliases; if FALSE, includes previous and alias symbols.
#'
#' @return A vector of unique gene aliases.
#'
#' @export
check_gene_aliases <- function(hgnc, approved = TRUE) {
  # Filter the mutmatch_geneAliases data frame to find any row where the hgnc symbol appears
  aliases <- mutmatch_geneAliases %>%
    dplyr::filter_all(dplyr::any_vars(stringr::str_detect(., base::paste0("^", hgnc, "$"))))

  # Further filter the aliases based on the specific columns where hgnc appears
  aliases <- aliases[`Approved symbol` == hgnc |
    `Previous symbols` == hgnc |
    hgnc %in% stringr::str_trim(unlist(stringr::str_split(`Alias symbols`, ",")))]

  # If only approved symbols are requested
  if (approved) {
    aliases <- aliases$`Approved symbol` %>%
      base::unique() %>%
      base::setdiff("")
  } else {
    # If all types of aliases are requested
    aliases <- c(
      aliases$`Approved symbol`,
      aliases$`Previous symbols`,
      stringr::str_trim(unlist(stringr::str_split(aliases$`Alias symbols`, ",")))
    ) %>%
      base::unique() %>%
      base::setdiff("")
  }

  return(aliases)
}
