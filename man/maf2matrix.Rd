\name{maf2matrix}
\alias{maf2matrix}

\title{
Convert a MAF file to a mutation matrix
}
\description{
Create a mutation matrix using variant types
}
\usage{
maf2matrix(maf, surv.data = NULL, sample.name = "Tumor_Sample_Barcode",
           gene.name = "Hugo_Symbol", variant.type = "Variant_Classification")
}

\arguments{
  \item{maf}{
a MAF file
}
  \item{surv.data}{
survival data for sample names (default: NULL)
}
  \item{sample.name}{
the column name containing sample names (defult: "Tumor_Sample_Barcode")
}
  \item{gene.name}{
the column name containing gene names  (defult: "Hugo_Symbol")
}
  \item{variant.type}{
the column name containing variant types  (defult: "Variant_Classification")
}
}


\value{
a mutation matrix
}


\author{
Dongmin Jung
}



\examples{
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml.maf <- read.csv(laml.maf, sep = "\t")
laml.mat <- maf2matrix(laml.maf)
}









