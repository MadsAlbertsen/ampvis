\name{amp_heatmap}
\alias{amp_heatmap}
\title{Generates a ggplot2 style heatmap from amplicon data in phyloseq format}
\usage{
amp_headtmap(data)
}
\arguments{
  \item{data}{(required) A phyloseq object including sample
  data.}

  \item{group}{A variable from the associated sample data
  to group samples by.}

  \item{scale}{A variable from the associated sample data
  to scale the abundance by.}

  \item{normalise}{A sample or group to normalise the
  counts to.}

  \item{tax.aggregate}{The taxonomic level that the data
  should be aggregated to (defualt: Phylum)}

  \item{tax.show}{The number of taxa to show or a vector of
  taxa names (default: 10).}

  \item{tax.clean}{Replace the phylum Proteobacteria with
  the respective Classes instead (default: T).}

  \item{order.x}{A taxonomy group or vector to order the
  x-axis by.}

  \item{order.y}{A sample or vector to order the y-axis
  by.}

  \item{plot.numbers}{Plot the values on the heatmap
  (default: T)}

  \item{plot.breaks}{A vector of breaks for the abundance
  legend.}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000)}

  \item{output}{To output a plot or the complete data
  inclusive dataframes (default: plot)}
}
\value{
A ggplot2 object or a list with the plot and associated
dataframes.
}
\description{
A nice long description
}
\author{
Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

