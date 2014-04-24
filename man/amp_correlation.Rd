\name{amp_correlation}
\alias{amp_correlation}
\title{Generates correlation networks from amplicon data.}
\usage{
amp_correlation(data)
}
\arguments{
  \item{data}{(required) A phyloseq object including sample
  data.}

  \item{trans}{Transform the raw counts, currently only
  supports "none" or "sqrt" (default: "none").}

  \item{treshold.cor}{Absolute correlation treshold
  (default: 0.8).}

  \item{treshold.count}{Minimum count in order to use a
  pair of observations for test of correlation (default:
  10)}

  \item{treshold.pval}{P-value treshold for correlations
  (default: 0.01).}

  \item{treshold.abundance}{Average abundance treshold
  (default: 0)}

  \item{ignore.zero}{Remove observations where 1
  obeservation has 0 counts from the correlation test
  (default: F)}

  \item{tax.clean}{Assign best classification to OTUs
  (default: T).}

  \item{label}{Plot taxonomic classifications (e.g.
  "Genus") instead of points (default: "OTU")}

  \item{scale.abundance}{Scale the size of nodes by
  abundance (default: F).}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000)}

  \item{scale.size}{Scale the size of the plotted objects
  (default: 1)}

  \item{output}{Either plot or complete (default: "plot").}
}
\value{
A ggplot2 object
}
\description{
A nice long description
}
\author{
Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

