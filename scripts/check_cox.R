# Python Cox fit vs survival::coxph, PH test
library(survival)
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1], check.names = FALSE)
x <- log(1e-3 + d[["CD4+PD-1+|CN-9"]])
f <- coxph(Surv(os, death) ~ x, data = d, ties = "breslow")
s <- summary(f)
cat(sprintf("beta %.4f  se %.4f  wald p %.4f  LR p %.4f\n", coef(f), sqrt(vcov(f)), s$coefficients[5], s$logtest["pvalue"]))
print(cox.zph(f))
