# qPCR analysis script

## load required packages
library(dplyr)

## import data
data <- read.csv("sampleAnalysis.csv")

## Test for normal distribution using Shapiro-Wilk's test
shapiro.test(data$MNE)

### assume data is normally distributed if p > 0.05. 
### normal distribution, use parametric tests
### non-normal distribution, use non-parametric test

## Paired t-test with Welch's correction
t.test(MNE ~ Group, data)