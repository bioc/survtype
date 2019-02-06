\name{quantile_normalization}
\alias{quantile_normalization}

\title{
Normalize a gene expression profile
}
\description{
Normalize expression data using quantile normalization
}
\usage{
quantile_normalization(x)
}

\arguments{
  \item{x}{
an expression profile
}
}


\value{
a normalized matrix
}


\author{
Dongmin Jung
}



\examples{
set.seed(1)
nrows <- 50
ncols <- 5
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
normalized.matrix <- quantile_normalization(counts)
}
