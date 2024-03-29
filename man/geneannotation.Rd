% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation.R
\name{geneannotation}
\alias{geneannotation}
\title{Annotate genes}
\usage{
geneannotation(genesymbols, geneentrezs)
}
\arguments{
\item{genesymbols}{A vector containing the symbols of the genes need to be
annotated.}

\item{geneentrezs}{A vector containing the ENTREZ IDs of the genes need to
be annotated. If \code{genesymbols} is NULL, this value will be used, and
if this value is also NULL, the function will turn NULL.}
}
\value{
A data.frame with the gene annotation information, including the
genomic coordinates of the probes and their function descriptions.
}
\description{
Annotate genes coordinates and functions
}
