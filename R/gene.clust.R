
gene.clust <- function(object, K, ...)
{
  Group <- object$cluster
  annotation.col <- data.frame(as.factor(Group))
  colnames(annotation.col) <- "Group"
  rownames(annotation.col) <- colnames(object$exprs.data)
  annotation.colors <- list(Group = c(scales::hue_pal()(K)))
  names(annotation.colors$Group) <- seq_len(K)
  # assign clusters for genes
  Cluster <- cutree(object$heatmap$tree_row, K)
  for(i in seq_len(K))
  {
    pheatmap(object$exprs.data[Cluster == i, object$heatmap$tree_col$order],
             cluster_cols = FALSE, annotation_col = annotation.col,
             annotation_colors = annotation.colors, ...)
  }
  Cluster
}
