\name{Single.survgroup}
\alias{Single.survgroup}

\title{
Patient group identification via survival information and expression of a single gene
}
\description{
All midpoints of the expression level or real-valued statistic are investigated to find the best threshold giving the most significant difference between two groups. Any patients having the value greater than the best cutoff are assigned as the "high-score" class. Otherwise, the others belong to the "low-score" class.
}
\usage{
Single.survgroup(surv.data, time, status, single.gene, intermediate = FALSE,
                 group.names = c("High", "Intermediate", "Low"))
}

\arguments{
  \item{surv.data}{
survival data
}
  \item{time}{
survival time
}
  \item{status}{
status indicator
}
  \item{single.gene}{
expression level of a single gene
}
  \item{intermediate}{
a logical value indicating whether or not the intermediate class is considered (default: FALSE)
}
  \item{group.names}{
the name of clusters for "high-score", "intermediate-score", and "low-score" classes (defult: "High", "Intermediate", "Low")
}
}



\value{
  \item{time}{
survival time
}
  \item{status}{
status indicator
}
  \item{surv.data}{
survival data
}
  \item{score}{
a vector of scores
}
  \item{summary}{
a list of thresholds, chi-squared statistics and p-values
}
  \item{cluster}{
a vector of clusters to which samples are assigned
}
  \item{fit}{
fitted survival curves
}
}



\author{
Dongmin Jung
}



\seealso{
survival::Surv, survival::survfit, survival::survdiff
}

\examples{
set.seed(1)
nrows <- 1
ncols <- nrow(ovarian)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colnames(counts) <- paste("X", 1:ncols, sep = "")
rownames(ovarian) <- paste("X", 1:ncols, sep = "")
Single.ovarian <- Single.survgroup(ovarian, time = "futime", status = "fustat", counts[1,])
plot(Single.ovarian, pval = TRUE)
}

