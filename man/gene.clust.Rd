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
DLBCLgenes <- read.csv("https://doi.org/10.1371/journal.pbio.0020108.sd012", header = FALSE)
DLBCLpatients <- read.csv("https://doi.org/10.1371/journal.pbio.0020108.sd013", header = FALSE)
colnames(DLBCLpatients) <- c("time", "status")
rownames(DLBCLpatients) <- colnames(DLBCLgenes)
DLBCL.survtype <- Exprs.survtype(DLBCLpatients, time = "time", status = "status",
                                 DLBCLgenes, num.genes = 50, scale = "row",
                                 gene.sel = TRUE, clustering_method = "ward.D2")
plot(DLBCL.survtype, pval = TRUE)
gene.clust(DLBCL.survtype, 2, scale = "row", clustering_method = "ward.D2")
}

