% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_annotation.R
\name{load_annotation}
\alias{load_annotation}
\title{Load the file(s) with annotation}
\usage{
load_annotation(
  annotationGenomeWidePath,
  annotationGenePath = NULL,
  hgnc = NULL
)
}
\arguments{
\item{annotationGenomeWidePath}{A path to the file with genome-wide annotation for each sample}

\item{annotationGenePath}{A path to the file with gene-specific annotation for each sample.
Can be 'sqlite' or 'db' if the file is produced with make_sql_gene_annotation function (recommended for large files). Default is NULL.}

\item{hgnc}{The HGNC symbol of the target gene in case a gene-specific annotation is used. Default is NULL.}
}
\value{
A data.table object
}
\description{
Load the file(s) with annotation
}
\examples{
\dontrun{
annotationGenomeWidePath <- system.file("extdata",
  "example_genomewide_annotation.csv",
  package = "mutmatch"
)
annotationGenePath <- system.file("extdata", "example_gene_annotation.csv.gz",
  package = "mutmatch"
)
hgnc <- "KRAS"
}
}
