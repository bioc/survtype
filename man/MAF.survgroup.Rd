\name{MAF.survgroup}
\alias{MAF.survgroup}

\title{
Patient group identification via survival data and mutation annotation information
}
\description{
The groups of patients are identified according to whether the varints exist on a single gene. Survival difference between patients with and without mutations is compared.
}
\usage{
MAF.survgroup(surv.data, time, status, maf, variants = NULL,
              sample.name = "Tumor_Sample_Barcode", gene.name = "Hugo_Symbol",
              variant.type = "Variant_Classification", num.genes = 10,
              siginificant.genes = 1, ...)
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
  \item{maf}{
a MAF file
}
  \item{variants}{
types of varints on a single gene for mutated samples. samples with any mutations, defined as mutated samples by default.
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
  \item{num.genes}{
the number of top genes based on the number of mutated genes (default: 10)
}
  \item{siginificant.genes}{
the number of top genes based on the statistical siginificance of mutated genes (default: 1)
}
  \item{...}{
additional parameters for the "ggsurvplot" for the statistically siginificant genes
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
  \item{maf.matrix}{
a mutation matrix
}
  \item{summary}{
a list of number of samples with variants, chi-squared statistics and p-values
}
  \item{cluster}{
a vector of integers indicating the cluster to which each sample is assigned, for the most significant gene
}
  \item{fit}{
fitted survival curves, for the most significant gene
}
}

\author{
Dongmin Jung
}



\seealso{
survival::Surv, survival::survfit, survival::survdiff, survminer::ggsurvplot
}

\examples{
laml.maf <- system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
laml.maf <- read.csv(laml.maf, sep = "\t")
laml.clinical.data <- read.csv(laml.clin, sep = "\t", row.names = 1)
index <- which(laml.clinical.data$days_to_last_followup == -Inf)
laml.clinical.data <- laml.clinical.data[-index,]
laml.clinical.data <- data.frame(laml.clinical.data)
laml.survgroup <- MAF.survgroup(laml.clinical.data, time = "days_to_last_followup",
                                status = "Overall_Survival_Status", laml.maf,
                                num.genes = 5, siginificant.genes = 2,
                                pval = TRUE)
}

