#' Make SQL Gene Annotation Database
#'
#' This function reads gene annotation data from a specified file and writes it to an SQLite database.
#'
#' @param annotationGenePath A character string specifying the path to the file containing gene annotation data.
#' The file at this path will be read into R.
#'
#' @param annotationGenePathSQL A character string specifying the path where the SQL database will be stored.
#' This path must end with a '.sqlite' or '.db' extension to indicate that it is an SQLite database file.
#'
#' @examples
#' \dontrun{
#' make_sql_gene_annotation(
#'   annotationGenePath = system.file("extdata",
#'     "example_gene_annotation.csv.gz",
#'     package = "mutmatch"
#'   ),
#'   annotationGenePathSQL = "/path/to/example_gene_annotation.sqlite"
#' )
#' }
#' @export
make_sql_gene_annotation <- function(annotationGenePath, annotationGenePathSQL) {
  # Check the file extension and modify if needed
  if (!grepl("\\.(sqlite|db)$", annotationGenePathSQL)) {
    stop("Invalid file extension for annotationGenePathSQL. Must be '.sqlite' or '.db'.")
  }

  # Read the data
  sample_annotation_gene_data <- data.table::fread(annotationGenePath)

  # Connect to a new SQLite database at the specified path
  # (it will be created if it doesn't exist)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = annotationGenePathSQL)

  # Write the data.frame to the SQLite database
  DBI::dbWriteTable(con, "sample_annotation_gene_data", sample_annotation_gene_data, overwrite = TRUE)

  # Create an index on the Gene column to speed up queries
  DBI::dbExecute(con, "CREATE INDEX idx_gene ON sample_annotation_gene_data(Gene)")

  # Disconnect from the SQLite database
  DBI::dbDisconnect(con)
}
