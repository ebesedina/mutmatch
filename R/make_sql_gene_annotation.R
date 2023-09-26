make_sql_gene_annotation = function(annotationGenePath) {
  # Read the data
  sample_annotation_gene_data <- data.table::fread(annotationGenePath)

  # Connect to a new SQLite database (it will be created if it doesn't exist)
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "gene_database.sqlite")

  # Write the data.frame to the SQLite database
  DBI::dbWriteTable(con, "sample_annotation_gene_data", sample_annotation_gene_data, overwrite = TRUE)

  # Disconnect from the SQLite database
  DBI::dbDisconnect(con)
}

