
plot <- function(object, ...) UseMethod("plot")

plot.survtype <- function(object, ...)
{
  Time <- object$time
  Status <- object$status
  Group <- object$cluster
  surv.data <- data.frame(Time, Status, Group)
  fit <- object$fit
  if(!is.null(object$subset))
  {
    Time <- Time[object$subset]
    Status <- Status[object$subset]
    Group <- Group[object$subset]
    surv.data <- data.frame(Time, Status, Group)
  }
  ggsurvplot(fit, data = surv.data, ...)
}
