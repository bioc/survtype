
MAF.survgroup <- function(surv.data, time, status, maf, variants = NULL,
                          sample.name = "Tumor_Sample_Barcode", 
                          gene.name = "Hugo_Symbol",
                          variant.type = "Variant_Classification",
                          num.genes = 10, siginificant.genes = 1, ...)
{
  maf.matrix <- maf2matrix(maf, surv.data, sample.name = sample.name,
                           gene.name = gene.name, variant.type = variant.type)
  surv.data <- surv.data[colnames(maf.matrix),]
  if (!identical(rownames(surv.data), colnames(maf.matrix)))
    stop("the sample ID of survivla data and MAF is not identical")
  if(all(variants %in% unique(maf[,variant.type])) == FALSE)
    stop("some of variant types to be considered are not valid")
  Time <- surv.data[,time]
  Status <- surv.data[,status]
  Status.levels <- unique(Status)
  Status[Status == Status.levels[1]] <- 0
  Status[Status == Status.levels[2]] <- 1
  # find p-values for selected genes
  num.samples.with.variants <- 0
  p.values <- 0
  chisq.stat <- 0
  pattern <- ifelse(is.null(variants),
                    paste(unique(maf[,variant.type]), collapse = "|"),
                    paste(variants, collapse = "|"))
  for(i in seq_len(num.genes))
  {
    index <- grepl(pattern, maf.matrix[i,])
    num.samples.with.variants[i] <- sum(index)
    Group <- as.numeric(index)
    chisq.stat[i] <- survdiff(Surv(Time, Status) ~ Group,
                              cbind(surv.data, Group))$chisq
    p.values[i] <- ifelse(sum(index) == 0, NA,
                          1 - pchisq(chisq.stat[i], 1))
  }
  output <- NULL
  output$time <- Time
  output$status <- Status
  output$surv.data <- surv.data
  output$maf.matrix <- maf.matrix
  output$summary <- cbind(num.samples.with.variants, chisq.stat, p.values)
  rownames(output$summary) <- rownames(maf.matrix)[seq_len(num.genes)]
  output$summary <- data.frame(output$summary[order(p.values),])
  for(i in seq_len(siginificant.genes))
  {
    Group <- NULL
    Group <- ifelse(as.matrix(grepl(pattern, maf.matrix[rownames(output$summary)[i],])) == 1,
                    paste(rownames(output$summary)[i], "MT"),
                    paste(rownames(output$summary)[i], "WT"))
    fit <- survfit(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
    if(i == 1)
    {
      output$cluster <- Group[,1]
      output$fit <- fit
    }
    p <- ggsurvplot(fit, data = data.frame(Time, Status, Group), ...)
    print(p)
  }
  class(output) <- c("survtype")
  output
}
