#' Get Annotated Mutations
#'
#' This function identifies mutations from grMutations that are located inside grObject, fetches nucleotide sequences,
#' and creates contextual mutation strings.
#'
#' @param grObject Genomic range object specifying the regions of interest.
#' @param grMutations Genomic range object containing mutation data.
#' @param nt_length A character string specifying the nucleotide length for context.
#'                  Can be either "tri" or "penta".
#'
#' @return A data.table containing the annotated mutations. If there are no overlaps,
#'         returns NULL.
#'
#' @examples
#' # grObject and grMutations should be GenomicRanges objects
#' # Here is how you would typically call the function:
#' # result <- get_mutations_annot(grObject, grMutations, nt_length = "tri")
#'
#' @export
get_mutations_annot <- function(grObject, grMutations, nt_length = c("tri", "penta")) {
  # Set the `nt_length` to a validated value ('tri' or 'penta') and determine the offset value
  nt_length <- base::match.arg(nt_length)
  offset <- base::ifelse(nt_length == "tri", 1, 2)

  # Ensure the genomic ranges are using the UCSC style
  GenomeInfoDb::seqlevelsStyle(grMutations) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(grObject) <- "UCSC"

  # Find overlaps between mutations and the object, ignoring strand info
  mutationsData <- GenomicRanges::findOverlaps(query = grMutations, subject = grObject, type = "within", ignore.strand = TRUE)

  # If there are overlapping mutations
  if (base::length(mutationsData) > 0) {
    # Subset the overlapping mutations and convert to a data.table
    mutationsData <- grMutations[S4Vectors::queryHits(mutationsData), ] %>%
      data.table::as.data.table() %>%
      dplyr::select(-c("strand", "width"))

    # Remove any duplicate mutations
    mutationsData <- mutationsData %>% base::unique()

    # Fetch the surrounding sequence context for each mutation
    mutationsData[, Context := {
      seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, seqnames, start - offset, end + offset)
      base::as.character(seqs)
    }]

    # Adjust the sequence context for germline mutations
    mutationsData[, Context := base::paste0(substr(Context, 1, offset), ref_allele, substr(Context, offset + 2, offset * 2 + 1))]
    mutationsData[, Mutation := base::paste(Context, mutated_to_allele, sep = ">")]

    # For mutations of A or G, reverse-complement the sequence context
    rev_complement <- mutationsData[ref_allele %in% c("A", "G")]

    # If there are any such mutations
    if (nrow(rev_complement) > 0) {
      # Get reverse complement of the sequence context
      rev_complement[, Context := as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(Context)))]
      rev_complement[, Mutation := paste(Context, chartr("ATGC", "TACG", mutated_to_allele), sep = ">")]

      # Update these in the main mutationsData
      mutationsData[ref_allele %in% c("A", "G")] <- rev_complement
    }

    # Return the mutations without the "Context" column
    return(mutationsData[, -"Context"])
  }

  # If there are no overlapping mutations, return an empty data.table with the correct column names
  mutationsDataCols <- grMutations %>%
    data.table::as.data.table() %>%
    dplyr::select(-c("strand", "width")) %>%
    base::colnames() %>%
    c("Mutation")

  mutationsData <- data.table::data.table(base::matrix(ncol = length(mutationsDataCols), nrow = 0))
  data.table::setnames(mutationsData, mutationsDataCols)
  return(mutationsData)
}
