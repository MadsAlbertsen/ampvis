\name{amp_point}
\alias{amp_point}
\title{Generates a ggplot2 style plot of the most abundant OTUs}
\usage{
amp_point(data)
}
\arguments{
  \item{data}{(required) A phyloseq object including sample
  data.}

  \item{tax.show}{A number of taxa to show or a vector of
  taxa names (default: 10).}

  \item{tax.clean}{Replace the phylum Proteobacteria with
  the respective Classes instead (default: T).}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000)}

  \item{plot.type}{Either point or boxplot
  (default:point).}
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

