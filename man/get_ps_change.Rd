% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Downstream Analysis.R
\name{get_ps_change}
\alias{get_ps_change}
\title{Create a data frame with the probability of differential expression for each pause site}
\usage{
get_ps_change(sites)
}
\arguments{
\item{sites}{List of ribosome pause sites; usually the output of \code{get_pause_sites()}}
}
\value{
Data frame with one row for each pause site and columns containing a unique identifier,
transcript location, nucleotide position, log fold change (logFC), and p-value of change across conditions for each pause site.
}
\description{
Takes as input a list of ribosome pause sites (the output of \code{get_pause_sites()}) and returns a
data frame with one row for each pause site and columns containing information about pause site location
and differential expression (p-value, logFC).
}
\details{
The function \code{\link{get_ps_change}} calls edgeR's methods to quantify change across conditions
with a limited number of samples, producing a log fold-change (logFC) and p-value. You can read more about these methods
\href{https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf}{here}.
}
