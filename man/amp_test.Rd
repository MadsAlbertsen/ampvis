\name{amp_test}
\alias{amp_test}
\title{Tests if there is a significant difference in abundance between selected conditions}
\usage{
amp_test(data, group)
}
\arguments{
  \item{data}{(required) A phyloseq object including sample
  data.}

  \item{group}{(required) A variable from the associated
  sample data to group samples by.}

  \item{tax.aggregate}{Group data at specific taxonomic
  level (defaul: "Genus").}

  \item{tax.clean}{Replace the phylum Proteobacteria with
  the respective Classes instead (default: T).}

  \item{tax.empty}{Either "remove" OTUs without taxonomic
  information or "rename" with OTU ID (default: rename).}

  \item{scale.seq}{The number of sequences in the
  pre-filtered samples (default: 20000).}

  \item{output}{To output the "pval" or the "complete" data
  inclusive dataframes (default: "pval").}

  \item{sig}{Significance cutoff to report results
  (default: 1).}
}
\value{
A p-value for each comparison.
}
\description{
A nice long description
}
\author{
Mads Albertsen \email{MadsAlbertsen85@gmail.com}
}

