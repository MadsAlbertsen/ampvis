\name{amp_ordinate}
\alias{amp_ordinate}
\title{Generates a ggplot2 style ordinate plot of amplicon data.}
\usage{
amp_ordinate(data)
}
\arguments{
  \item{data}{(required) A phyloseq object including sample
  data.}

  \item{trans}{Transform the raw counts, currently only
  supports "none" or "sqrt" (default: "sqrt").}

  \item{ordinate.type}{Either PCA or NMDS (default:
  "PCA").}

  \item{ncomp}{The number of principal components to
  extract (default: 5)}

  \item{plot.x}{Variable to plot on the x-axis (default:
  "PC1")}

  \item{plot.y}{Variable to plot on the y-axis (default:
  "PC2")}

  \item{plot.color}{Color the points by a sample variable.}

  \item{plot.point.size}{Size of the plotted points
  (default: 3)}

  \item{plot.species}{Plot loadings as points (default: F)}

  \item{plot.nspecies}{Plot the n most extreme species with
  their genus classification (default: 0).}

  \item{envfit.factor}{A vector of factor variables from
  the sample data used for envfit to the model.}

  \item{envfit.numeric}{A vector of numerical variables
  from the sample data used for envfit to the model.}

  \item{envfit.significant}{The significance treshold for
  displaying envfit parameters (default: 0.01).}

  \item{envfit.resize}{Scale the size of the numeric arrows
  (default: 1).}

  \item{tax.clean}{Add best assignment to Genus level
  classification if none exists.}

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

