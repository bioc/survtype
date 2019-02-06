\name{gene.clust}
\alias{gene.clust}

\title{
Plots of the heatmap for each cluster of expression data
}
\description{
Heatmaps of clustered genes for subtypes of samples can be drawn.
}
\usage{
gene.clust(object, K, ...)
}

\arguments{
  \item{object}{
the result of "Exprs.survtype"
}
  \item{K}{
the number of clusters
}
  \item{...}{
additional parameters for the "pheatmap"
}
}


\value{
Heatmap for each cluster
}


\author{
Dongmin Jung
}



\seealso{
pheatmap::pheatmap
}

\examples{
set.seed(1)
nrows <- 200
ncols <- nrow(ovarian)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colnames(counts) <- paste("X", 1:ncols, sep = "")
rownames(ovarian) <- paste("X", 1:ncols, sep = "")
SE <- SummarizedExperiment(assays = SimpleList(counts = counts))
ovarian.survtype <- Exprs.survtype(ovarian, time = "futime", status = "fustat",
                                 assay(SE), num.genes = 50, scale = "row",
                                 clustering_method = "ward.D2")
plot(ovarian.survtype, pval = TRUE)
gene.clust(ovarian.survtype, 2, scale = "row", clustering_method = "ward.D2")
}

