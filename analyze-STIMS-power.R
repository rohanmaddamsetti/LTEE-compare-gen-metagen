## analyze-STIMS-power.R by Rohan Maddamsetti and Nkrumah Grant. 
## This script takes the output of generate-STIMS-power-data.jl,
## and makes figures.

## References:
## https://en.wikipedia.org/wiki/Sensitivity_and_specificity
## https://en.wikipedia.org/wiki/Power_of_a_test
## https://en.wikipedia.org/wiki/Confusion_matrix

library(tidyverse)
library(cowplot)


## Count the number of True Positives and False Negatives when testing
## the positive selection module for positive selection, and
## count the number of False Positives and True Negatives when testing
## the neutral module for positive selection.
make.positive.confusion.matrix <- function(power.analysis.df, pval.threshold=0.05) {
    ## assume positive selection module here, testing for positive selection.
    positive.part <- power.analysis.df %>%
        filter(GeneModule == "positive") %>%
    group_by(MutatorStatus, Variable, VariableValue) %>%
        mutate(Outcome = ifelse(Pval < pval.threshold, T, F)) %>%
        summarize(TP = sum(Outcome == T), FN = n() - TP)
    
    ## assume neutral module here, testing for positive selection.
    neutral.part <- power.analysis.df %>%
        filter(GeneModule == "neutral") %>%
    group_by(MutatorStatus, Variable, VariableValue) %>%
        mutate(Outcome = ifelse(Pval < pval.threshold, T, F)) %>%
        summarize(FP = sum(Outcome == T), TN = n() - FP)
    
    ## now, merge these two dataframes to make the confusion matrix.
    confusion.matrix <- full_join(positive.part, neutral.part) %>%
        ## See the definitions here:
        ## https://en.wikipedia.org/wiki/Sensitivity_and_specificity
        ## Statistical power is a synonym for sensitivity.
        mutate(Sensitivity = 1 - (FN / (TP + FN))) %>%
        mutate(Specificity = 1 - (FP / (FP + TN)))
    return(confusion.matrix)
}


## Count the number of True Positives and False Negatives when testing
## the purifying selection module for purifying selection, and
## count the number of False Positives and True Negatives when testing
## the neutral module for purifying selection.
make.negative.confusion.matrix <- function(power.analysis.df, pval.threshold=0.05) {
    ## assume negative selection module here, testing for negative selection.
    negative.part <- power.analysis.df %>%
        filter(GeneModule == "purifying") %>%
    group_by(MutatorStatus, Variable, VariableValue) %>%
        mutate(Outcome = ifelse(Pval > (1 - pval.threshold), T, F)) %>%
        summarize(TP = sum(Outcome == T), FN = n() - TP)
    
    ## assume neutral module here, testing for positive selection.
    neutral.part <- power.analysis.df %>%
        filter(GeneModule == "neutral") %>%
    group_by(MutatorStatus, Variable, VariableValue) %>%
        mutate(Outcome = ifelse(Pval > (1 - pval.threshold), T, F)) %>%
        summarize(FP = sum(Outcome == T), TN = n() - FP)
    
    ## now, merge these two dataframes to make the confusion matrix.
    confusion.matrix <- full_join(negative.part, neutral.part) %>%
        ## See the definitions here:
        ## https://en.wikipedia.org/wiki/Sensitivity_and_specificity
        ## Statistical power is a synonym for sensitivity.
        mutate(Sensitivity = 1 - (FN / (TP + FN))) %>%
        mutate(Specificity = 1 - (FP / (FP + TN)))
    return(confusion.matrix)
}


## get the power analysis data.
power.analysis.df <- read.csv("../results/SLiM-results/STIMS-power-analysis-data.csv") %>%
    mutate(GeneModule = factor(GeneModule, levels = c("positive", "neutral", "purifying")))

positive.confusion.matrix <- make.positive.confusion.matrix(power.analysis.df)

negative.confusion.matrix <- make.negative.confusion.matrix(power.analysis.df)
######################################################################################
