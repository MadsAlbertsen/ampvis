\name{amp_rabund}
\alias{amp_rabund}
\title{Generates a ggplot2 style rank abundance plot of the most abundant OTUs}
\usage{
amp_rabund(data)
}
\arguments{
  \item{data}{(required) A phyloseq object.}

  \item{group}{Group the data based on a sample variable.}

  \item{tax.show}{The number of taxa to show or a vector of
  taxa names (default: 10).}

  \item{tax.clean}{Replace the phylum Proteobacteria with
  the respective Classes instead (default: T).}

  \item{tax.aggregate}{The taxonomic level that the data
  should be aggregated to (defualt: OTU)}

  \item{tax.empty}{Either "remove" OTUs without taxonomic
  information or "rename" with OTU ID (default: rename).}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000)}

  \item{plot.type}{Either point or boxplot
  (default:point).}

  \item{name2}{A taxonomic level used for naming (default:
  "Phylum").}

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

