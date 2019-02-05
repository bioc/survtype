
Exprs.survtype <- function(surv.data, time, status, exprs.data, K = 2,
                           num.genes = 100, gene.sel = FALSE,
                           gene.sel.opt = list(verbose = FALSE), ...)
{
  if(!identical(rownames(na.omit(surv.data)), colnames(exprs.data)))
  {
    cat("the row names of the survivla data and the column names of the expression data are not identical \n")
    matched.samples <- intersect(rownames(na.omit(surv.data)), colnames(exprs.data))
    if(length(matched.samples) == 0)
    {
      stop("no matched samples beween survival data and expression data")
    }
    surv.data <- surv.data[matched.samples,]
    exprs.data <- exprs.data[,matched.samples]
  }
  cat(paste(nrow(surv.data)), "matched samples")
  
  Time <- surv.data[,time]
  Status <- surv.data[,status]
  Status.levels <- unique(Status)
  Status[Status == Status.levels[1]] <- 0
  Status[Status == Status.levels[2]] <- 1
  # scoring variables
  scores <- apply(exprs.data, 1, function(x) coxph(Surv(Time, Status) ~ as.numeric(x))$score)
  names(scores) <- rownames(exprs.data)
  variables <- names(head(sort(scores, decreasing = TRUE), min(num.genes, length(scores))))
  # variable selection
  if(gene.sel == TRUE)
  {
    clust.var.sel <- do.call(clustvarsel, 
                             c(list(t(exprs.data[variables,]), K), gene.sel.opt))
    variables <- variables[clust.var.sel$subset]
  }
  pheatmap.silent <- pheatmap(exprs.data[variables,], silent = TRUE, ...)
  # assign clusters for samples
  Group <- cutree(pheatmap.silent$tree_col, K)
  fit <- survfit(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
  # heatmap
  annotation.col <- data.frame(as.factor(Group))
  colnames(annotation.col) <- "Group"
  rownames(annotation.col) <- colnames(exprs.data)
  annotation.colors <- list(Group = c(scales::hue_pal()(K)))
  names(annotation.colors$Group) <- seq_len(K)
  pheatmap(exprs.data[variables,], annotation_col = annotation.col,
           annotation_colors = annotation.colors, ...)
  res <- survdiff(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
  print(res)
  output <- NULL
  output$n <- res$n
  output$obs <- res$obs
  output$exp <- res$exp
  output$var <- res$var
  output$chisq <- res$chisq
  output$call <- res$call
  output$fit <- fit
  output$cluster <- Group
  names(output$cluster) <- rownames(surv.data)
  output$time <- Time
  output$status <- Status
  output$surv.data <- surv.data
  output$exprs.data <- exprs.data[variables,]
  output$heatmap <- pheatmap.silent
  class(output) <- c("survtype")
  output
}
