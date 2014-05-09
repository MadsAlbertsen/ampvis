\name{amp_core}
\alias{amp_core}
\title{Generates a ggplot2 style core community plots}
\usage{
amp_core(data)
}
\arguments{
  \item{data}{(required) A phyloseq object.}

  \item{group}{Group the data based on a sample variable
  (default: "Sample").}

  \item{tax.clean}{Replace the phylum Proteobacteria with
  the respective Classes instead (default: T).}

  \item{tax.aggregate}{The taxonomic level that the data
  should be aggregated to (defualt: OTU)}

  \item{tax.empty}{Either "remove" OTUs without taxonomic
  information or "rename" with OTU ID (default: rename).}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000)}

  \item{plot.type}{Either core or frequency (default:
  frequency).}

  \item{weight}{Weight the frequency by abundance (default:
  T).}

  \item{abund.treshold}{Treshold for considering something
  abundant in percent (default: 0.1).}

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

