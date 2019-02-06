\name{Exprs.survtype}
\alias{Exprs.survtype}

\title{
Sample subtype identification via survival information and gene expression data
}
\description{
For discovery of subtypes of samples that are both clinically relevant and biologically meaningful, the Cox regession model and hierarchical clustering are combined. 
}
\usage{
Exprs.survtype(surv.data, time, status, exprs.data, K = 2, num.genes = 100,
               gene.sel = FALSE, gene.sel.opt = list(verbose = FALSE), ...)
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
  \item{exprs.data}{
expression data
}
  \item{K}{
the number of clusters (default: 2)
}
  \item{num.genes}{
the number of top genes based on the Cox score, before variable selection (default: 100)
}
  \item{gene.sel}{
a logical value indicating whether or not gene selection for clustring is applied (default: FALSE)
}
  \item{gene.sel.opt}{
a list of options for the gene selection function "clustvarsel". "verbose" is set to FALSE as default.
}
  \item{...}{
additional parameters for the "pheatmap"
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
  \item{fit}{
fitted survival curves
}
  \item{cluster}{
a vector of integers indicating the cluster to which each sample is assigned
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
  \item{exprs.data}{
expression data
}
}


\references{
Bair, E., & Tibshirani, R. (2004). Semi-supervised methods to predict patient survival from gene expression data. PLoS biology, 2(4), e108.
}


\author{
Dongmin Jung
}



\seealso{
survival::Surv, survival::survfit, survival::survdiff, survival::coxph, clustvarsel::clustvarsel, pheatmap::pheatmap
}

\examples{
set.seed(1)
nrows <- 100
ncols <- nrow(ovarian)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colnames(counts) <- paste("X", 1:ncols, sep = "")
rownames(ovarian) <- paste("X", 1:ncols, sep = "")
SE <- SummarizedExperiment(assays = SimpleList(counts = counts))
ovarian.survtype <- Exprs.survtype(ovarian, time = "futime", status = "fustat",
                                 assay(SE), num.genes = 10, scale = "row",
                                 clustering_method = "ward.D2")
plot(ovarian.survtype, pval = TRUE)
}

