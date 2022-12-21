# qPCR analysis script

## load required packages
library(dplyr)
library(broom)


## import data
data <- read.csv("homozygotesAnalysis.csv")

## Test for normal distribution using Shapiro-Wilk's test
shapiro.test(data$MNE)

### assume data is normally distributed if p > 0.05. 
### normal distribution, use parametric tests
### non-normal distribution, use non-parametric test

## Paired t-test with Welch's correction
res <- tidy(t.test(MNE ~ Group, data))

## print t-test results
write.csv(res, "homozygotesResults.csv")

