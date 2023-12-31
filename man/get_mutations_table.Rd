% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutations_table.R
\name{get_mutations_table}
\alias{get_mutations_table}
\title{Retrieve mutations within a genomic range}
\usage{
get_mutations_table(grObject, grMutations, annotation = NULL)
}
\arguments{
\item{grObject}{Genomic range object specifying the regions of interest.}

\item{grMutations}{Genomic range object containing mutation data.}

\item{annotation}{Data table with annotation details for each sample.}
}
\value{
A data table containing the mutations within the genomic range,
stratified based on annotation and cluster number.
}
\description{
This function identifies mutations within a given genomic range (grObject)
and stratifies the result based on the annotation and cluster number.
}
\examples{
\dontrun{
# Assume grObject, grMutations, and annotation are properly defined
result <- get_mutations_table(grObject, grMutations, annotation)
}
}
