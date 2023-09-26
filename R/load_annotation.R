#' Load the file(s) with annotation
#'
#' @param annotationGenomeWidePath A path to the file with genome-wide annotation for each sample
#' @param annotationGenePath A path to the file with gene-specific annotation for each sample.
#' Can be 'sqlite' or 'db' if the file is produced with make_sql_gene_annotation function (recommended for large files). Default is NULL.
#' @param hgnc The HGNC symbol of the target gene in case a gene-specific annotation is used. Default is NULL.
#'
#' @return A data.table object
#' @export
#'
#' @examples
#' \dontrun{
#' annotationGenomeWidePath <- system.file("extdata",
#'   "example_genomewide_annotation.csv",
#'   package = "mutmatch"
#' )
#' annotationGenePath <- system.file("extdata", "example_gene_annotation.csv.gz",
#'   package = "mutmatch"
#' )
#' hgnc <- "KRAS"
#' }
load_annotation <- function(annotationGenomeWidePath,
                            annotationGenePath = NULL,
                            hgnc = NULL) {
  if (!base::is.null(annotationGenomeWidePath)) {
    sample_annotation_genome <- data.table::fread(annotationGenomeWidePath) %>% unique()
  }
  if (!base::is.null(annotationGenePath)) {
    # Retrieve possible secondary names for the gene of interest
    gene_names <- check_gene_aliases(hgnc = hgnc, approved = FALSE)

    # Check if file with gene annotation is sql DB file
    file_extension <- tools::file_ext(annotationGenePath)
    if (file_extension %in% c("sqlite", "db")) {
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = annotationGenePath)

      # Generate SQL pattern using LIKE and wildcards
      query_pattern <- paste(gene_names, collapse = "%' OR Gene LIKE '")
      query_pattern <- paste0("%", query_pattern, "%")

      # Write SQL query to fetch matching rows based on the column name where genes are stored
      query <- sprintf("SELECT * FROM sample_annotation_gene_data WHERE Gene LIKE '%s'", query_pattern)

      # Execute the query and fetch the results
      sample_annotation_gene <- DBI::dbGetQuery(con, query) %>% data.table::data.table()

      # Disconnect from the SQLite database
      DBI::dbDisconnect(con)

      sample_annotation_gene <- sample_annotation_gene[Gene %in% gene_names][, -"Gene"] %>% unique()
    } else {
      # Get only rows with gene name in it
      sample_annotation_gene <- data.table::fread(cmd = base::paste(
        "zcat",
        annotationGenePath,
        "| grep -E",
        base::paste0("'", base::paste(gene_names, collapse = "|"), "'")
      ))
      sample_annotation_gene_colnames <- base::scan(annotationGenePath, what = "", nlines = 1, sep = ",", quiet = T)
      colnames(sample_annotation_gene) <- sample_annotation_gene_colnames
      sample_annotation_gene <- sample_annotation_gene[Gene %in% gene_names][, -"Gene"] %>% unique()
    }
  }

  # Return or both annotations or only the one which is specified
  if (!base::is.null(annotationGenomeWidePath) &
    !base::is.null(annotationGenePath)) {
    sample_annotation <- base::merge(sample_annotation_genome, sample_annotation_gene, by = "Sample") %>% unique()
  } else if (!base::is.null(annotationGenomeWidePath) &
    base::is.null(annotationGenePath)) {
    sample_annotation <- sample_annotation_genome %>% unique()
  } else if (base::is.null(annotationGenomeWidePath) &
    !base::is.null(annotationGenePath)) {
    sample_annotation <- sample_annotation_gene %>% unique()
  }

  return(sample_annotation)
}
