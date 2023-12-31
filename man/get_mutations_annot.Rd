% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutations_annot.R
\name{get_mutations_annot}
\alias{get_mutations_annot}
\title{Get Annotated Mutations}
\usage{
get_mutations_annot(grObject, grMutations, nt_length = c("tri", "penta"))
}
\arguments{
\item{grObject}{Genomic range object specifying the regions of interest.}

\item{grMutations}{Genomic range object containing mutation data.}

\item{nt_length}{A character string specifying the nucleotide length for context.
Can be either "tri" or "penta".}
}
\value{
A data.table containing the annotated mutations. If there are no overlaps,
returns NULL.
}
\description{
This function identifies mutations from grMutations that are located inside grObject, fetches nucleotide sequences,
and creates contextual mutation strings.
}
\examples{
# grObject and grMutations should be GenomicRanges objects
# Here is how you would typically call the function:
# result <- get_mutations_annot(grObject, grMutations, nt_length = "tri")

}
