\name{Surv.survtype}
\alias{Surv.survtype}

\title{
Sample subtype identification via survival information
}
\description{
Any patient who lived longer than the median was considered to be a "low-risk" patient, and any patient that lived less than the median was considered to be a "high-risk"" patient. In this manner, we assigned a class label to each observation.
For censored data, we can estimate the probability that each censored observation belongs to the "low-risk" and "high-risk" classes, respectively.
}
\usage{
Surv.survtype(surv.data, time, status)
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
}


\value{
  \item{n}{
the number of subjects in each group
}
  \item{obs}{
the weighted observed number of events in each group
}
  \item{exp}{
the weighted expected number of events in each group
}
  \item{chisq}{
the chi-squared statistic for a test of equality
}
  \item{call}{
the matched call
}
  \item{cluster}{
a vector of clusters to which samples are assigned
}
  \item{time}{
survival time
}
  \item{status}{
status indicator
}
  \item{surv.data}{
survival data
}
  \item{fit}{
fitted survival curves
}
}


\references{
Bair, E., & Tibshirani, R. (2004). Semi-supervised methods to predict patient survival from gene expression data. PLoS biology, 2(4), e108.
}


\author{
Dongmin Jung
}



\seealso{
survival::Surv, survival::survfit, survival::survdiff
}

\examples{
data(ovarian)
ovarian.survtype <- Surv.survtype(ovarian, time = "futime", status = "fustat")
plot(ovarian.survtype, pval = TRUE)
}

