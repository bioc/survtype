
Surv.survtype <- function(surv.data, time, status)
{
  Median <- median(surv.data[,time])
  Time <- surv.data[,time]
  Status <- surv.data[,status]
  Status.levels <- unique(Status)
  Status[Status == Status.levels[1]] <- 0
  Status[Status == Status.levels[2]] <- 1
  fit <- survfit(Surv(Time, Status) ~ 1)
  high.risk.prob <- 1 - summary(fit, times = Time)$surv/
    summary(fit, times = Median)$surv
  # median-cut
  Group <- NA
  Group <- ifelse(Time >= Median, "Low Risk", Group)
  Group <- ifelse(is.na(Group) == TRUE & Status == 1 & Time < Median,
                  "High Risk", Group)
  Group <- ifelse(is.na(Group) == TRUE & Status == 0 &
                    Time < Median & high.risk.prob < 0.5,
                  "Low Risk", Group)
  Group <- ifelse(is.na(Group) == TRUE & Status == 0 &
                    Time < Median & high.risk.prob > 0.5,
                  "High Risk", Group)
  fit <- survfit(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
  res <- survdiff(Surv(Time, Status) ~ Group, cbind(surv.data, Group))
  print(res)
  output <- NULL
  output$n <- res$n
  output$obs <- res$obs
  output$exp <- res$exp
  output$var <- res$var
  output$chisq <- res$chisq
  output$call <- res$call
  output$cluster <- Group
  names(output$cluster) <- rownames(surv.data)
  output$time <- Time
  output$status <- Status
  output$surv.data <- surv.data
  output$fit <- fit
  class(output) <- c("survtype")
  output
}
