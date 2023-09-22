#' Load and Preprocess Mutation Data
#'
#' This function reads in mutation data from a specified file path and performs several preprocessing steps.
#' It offers the option to focus only on mutations for a specific chromosome.
#' The function also converts the mutation table to a GRanges object, retains important columns, and standardizes chromosome naming.
#' Finally, it filters out any SNVs that don't meet specific length criteria for the reference and mutated alleles.
#'
#' @param mutationsPath A path to the csv.gz file with mutations for each sample.
#' @param chr Optional. The chromosome for which mutations should be specifically loaded.
#' If not specified, mutations for all chromosomes are loaded.
#'
#' @return A GRanges object containing filtered and preprocessed mutations data.
#' @export
#'
#' @examples
#' \dontrun{
#' grSNVs <- load_mutations(mutationsPath = system.file("extdata",
#'   "example_mutations.csv.gz",
#'   package = "mutmatch"
#' ))
#'
#' grSNVs_chr1 <- load_mutations(mutationsPath = system.file("extdata",
#'   "example_mutations.csv.gz",
#'   package = "mutmatch"
#' ), chr = "1")
#' }
load_mutations <- function(mutationsPath,
                           chr = NULL) {
  # Load only mutations of a specific chromosome if it is specified
  if (!base::is.null(chr)) {
    snv_table <- data.table::fread(cmd = base::paste("zcat", mutationsPath, "| grep", base::paste0("'", chr, ",", "'")))
    snv_table_colnames <- base::scan(mutationsPath, what = "", nlines = 1, sep = ",", quiet = T)
    base::colnames(snv_table) <- snv_table_colnames
    snv_table <- snv_table[CHROM %in% c(chr, base::paste0("chr", chr))] %>% base::unique()
  } else {
    snv_table <- data.table::fread(mutationsPath)
    snv_table <- snv_table %>% base::unique()
  }

  # Convert to Granges object and leave only the important columns
  grSNVs <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(snv_table$CHROM),
    ranges = IRanges::IRanges(start = snv_table$POS, end = snv_table$POS),
    Sample = snv_table$Sample,
    ref_allele = snv_table$REF,
    mutated_to_allele = snv_table$ALT
  )

  # Standardize naming of the chromosomes
  GenomeInfoDb::seqlevelsStyle(grSNVs) <- "UCSC"

  # Filter SNVs
  grSNVs <- grSNVs[base::nchar(grSNVs$ref_allele) == 1 &
    base::nchar(grSNVs$mutated_to_allele) == 1]
  return(grSNVs)
}
