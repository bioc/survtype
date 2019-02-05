
Single.survgroup <- function(surv.data, time, status, single.gene, intermediate = FALSE,
                             group.names = c("High", "Intermediate", "Low"))
{
  if(!identical(rownames(na.omit(surv.data)), names(single.gene)))
  {
    cat("the row names of the survivla data and the names of the single gene are not identical \n")
    matched.samples <- intersect(rownames(na.omit(surv.data)), names(single.gene))
    if(length(matched.samples) == 0)
    {
      stop("no matched samples beween survival data and expression data")
    }
    surv.data <- surv.data[matched.samples,]
    single.gene <- single.gene[matched.samples]
  }
  cat(paste(nrow(surv.data)), "matched samples \n")
  
  Time <- surv.data[,time]
  Status <- surv.data[,status]
  Status.levels <- unique(Status)
  Status[Status == Status.levels[1]] <- 0
  Status[Status == Status.levels[2]] <- 1
  # midpoints
  single.gene <- as.numeric(single.gene)
  sort.single.gene <- sort(unique(single.gene))
  threshold <- approx(sort.single.gene, n = 2*length(sort.single.gene) - 1)$y
  threshold <- threshold [-which(threshold  %in% sort.single.gene)]
  # find p-values for all threshold
  if(intermediate == FALSE)
  {
    p.values <- 0
    chisq.stat <- 0
    for(i in seq_len(length(threshold)))
    {
      Group <- ifelse(single.gene > threshold[i], group.names[1], group.names[length(group.names)])
      SurvDiff <- survdiff(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
      chisq.stat[i] <- SurvDiff$chisq
      p.values[i] <- 1 - pchisq(SurvDiff$chisq, 1)
    }
    output <- NULL
    output$time <- Time
    output$status <- Status
    output$surv.data <- surv.data
    output$single.gene <- single.gene
    output$summary <- data.frame(cbind(threshold, chisq.stat, p.values)[order(p.values),])
    Group <- ifelse(single.gene > output$summary[1,1], group.names[1], group.names[length(group.names)])
    fit <- survfit(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
    output$cluster <- Group
    output$fit <- fit
  }
  if(intermediate == TRUE)
  {
    p.values <- matrix(NA, length(threshold), length(threshold))
    chisq.stat <- matrix(NA, length(threshold), length(threshold))
    num.Group1 <- matrix(NA, length(threshold), length(threshold))
    num.Group2 <- matrix(NA, length(threshold), length(threshold))
    num.Group3 <- matrix(NA, length(threshold), length(threshold))
    pb <- txtProgressBar(min = 0, max = length(threshold), style = 3)
    for(i in seq_len(length(threshold)))
    {
      for(j in seq_len(length(threshold)))
      {
        Group <- rep(group.names[2], length(single.gene))
        Group <- ifelse(single.gene > threshold[j], group.names[1], Group)
        Group <- ifelse(single.gene <= threshold[i], group.names[3], Group)
        num.Group1[i,j] <- sum(Group == group.names[1])
        num.Group2[i,j] <- sum(Group == group.names[2])
        num.Group3[i,j] <- sum(Group == group.names[3])
        Subset <- Group != group.names[2]
        SurvDiff <- survdiff(Surv(Time, Status) ~ Group, cbind(surv.data, Group),
                             subset = Subset)
        chisq.stat[i,j] <- SurvDiff$chisq
        p.values[i,j] <- 1 - pchisq(SurvDiff$chisq, 1)
      }
      Sys.sleep(0.01)
      setTxtProgressBar(pb, i)
    }
    cat("\n summarizing...")
    all.p.values <- unique(sort(p.values))
    index <- NULL
    for(i in seq_len(length(all.p.values)))
    {
      index <- rbind(index, which(p.values == all.p.values[i], arr.ind = TRUE))
    }
    
    output <- NULL
    output$time <- Time
    output$status <- Status
    output$surv.data <- surv.data
    output$summary <- data.frame(cbind(threshold[index[,2]], threshold[index[,1]],
                                       chisq.stat[index], p.values[index],
                                       num.Group1[index], num.Group2[index],
                                       num.Group3[index]))
    colnames(output$summary) <- c("threshold.high", "threshold.low", "chisq.stat",
                                  "p.values", "num.high", "num.intermediate", "num.low")
    Group <- rep(group.names[2], length(single.gene))
    Group <- ifelse(single.gene > output$summary[1,1], group.names[1], Group)
    Group <- ifelse(single.gene <= output$summary[1,2], group.names[3], Group)
    output$cluster <- Group
    Subset <- Group != group.names[2]
    output$subset <- Subset
    surv.data <- surv.data[Subset,]
    Time <- Time[Subset]
    Status <- Status[Subset]
    Group <- Group[Subset]
    fit <- survfit(Surv(Time, Status) ~ Group)
    output$fit <- fit
  }
  class(output) <- c("survtype")
  output
}
