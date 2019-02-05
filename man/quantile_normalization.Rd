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
DLBCLgenes <- read.csv("https://doi.org/10.1371/journal.pbio.0020108.sd012", header = FALSE)
normalized.matrix <- quantile_normalization(DLBCLgenes)
}
