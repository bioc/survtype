
maf2matrix <- function(maf, surv.data = NULL, sample.name = "Tumor_Sample_Barcode",
                       gene.name = "Hugo_Symbol", variant.type = "Variant_Classification")
{
  genes <- unique(maf[,gene.name])
  if(is.null(surv.data))
  {
    samples <- unique(maf[,sample.name])
  }
  else
  {
    samples <- intersect(unique(maf[,sample.name]),
                         rownames(surv.data))
  }
  # mutation matrix
  mutation <- matrix("", length(genes), length(samples))
  rownames(mutation) <- genes
  colnames(mutation) <- samples
  for(i in seq_len(ncol(mutation)))
  {
    index <- which(samples[i] == maf[,sample.name])
    sample.mut <- maf[index, c(gene.name, variant.type)]
    for(j in seq_len(nrow(sample.mut)))
    {
      index2 <- which(sample.mut[j, 1] == genes)
      mutation[index2, i] <- paste(mutation[index2, i], sample.mut[j, 2],
                                   ";", sep = "")
    }
  }
  mutation[order(rowSums(mutation != ""), decreasing = TRUE),]
}
