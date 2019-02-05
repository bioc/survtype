\name{plot.survtype}
\alias{plot.survtype}
\alias{plot}

\title{
Plot of survival curves of sample subtypes
}
\description{
Survival curves for subtypes of samples can be drawn.
}
\usage{
\method{plot}{survtype}(object, ...)
}

\arguments{
  \item{object}{
object of class "survtype"
}
  \item{...}{
additional parameters for the "ggsurvplot"
}
}


\value{
Survival curves
}


\author{
Dongmin Jung
}



\seealso{
survminer::ggsurvplot
}

\examples{
data(lung)
lung.survtype <- Surv.survtype(lung, time = "time", status = "status")
plot(lung.survtype, pval = TRUE)
}

