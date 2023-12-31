% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mutation_spectra.R
\name{get_mutation_spectra}
\alias{get_mutation_spectra}
\title{Generate Possible Mutation Spectra}
\usage{
get_mutation_spectra(nt_length = c("tri", "penta"))
}
\arguments{
\item{nt_length}{A character string specifying the nucleotide context length for constructing the mutation spectra.
Options are "tri" for trinucleotides (MS96) and "penta" for pentanucleotides. Default is c("tri", "penta").}
}
\value{
A data table with columns "Context" and "Mutation". The "Context" column contains the nucleotide
context in which the mutation occurs, and the "Mutation" column contains the mutation represented as "Context>MutatedNucleotide".
}
\description{
This function generates all possible mutation spectra given a nucleotide context length.
It constructs trinucleotide or pentanucleotide combinations, expands them with possible mutations,
and then filters the results to create a list of possible mutations.
}
\examples{
\dontrun{
get_mutation_spectra() # By default, uses trinucleotide context
get_mutation_spectra(nt_length = "penta") # Uses pentanucleotide context
}
}
