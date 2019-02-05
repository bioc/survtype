
quantile_normalization <- function(x)
{
  # x : expression profile
  Rownames <- rownames(x)
  Colnames <- colnames(x)
  y <- rowMeans(apply(x, 2, sort))
  qnorm <- apply(x, 2, function(x){quantile(y, ecdf(x)(x), type = 1)})
  dimnames(qnorm) <- list(Rownames, Colnames)
  qnorm
}
