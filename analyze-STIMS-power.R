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


## generic plotting function
.plot.sensitivity <- function(variable, xlabel, confusion.matrix, mutator.status) {
    confusion.matrix %>%
        filter(MutatorStatus == mutator.status) %>%
        filter(Variable == variable) %>%
        ggplot(aes(x = VariableValue, y = Sensitivity)) +
        theme_classic() +
        ylim(0,1) +
        xlab(xlabel) +
        geom_point() +
        geom_line()
}


## generic plotting function
.plot.specificity <- function(variable, xlabel, confusion.matrix, mutator.status) {
    confusion.matrix %>%
        filter(MutatorStatus == mutator.status) %>%
        filter(Variable == variable) %>%
        ggplot(aes(x = VariableValue, y = Specificity)) +
        theme_classic() +
        ylim(0,1) +
        xlab(xlabel) +
        geom_point() +
        geom_line()
}


## Specific plotting functions that are called.
plot.filtering_threshold.sensitivity <- partial(
    .f=.plot.sensitivity,
    "filtering_threshold",
    "Detection threshold (allele frequency)")


plot.filtering_threshold.specificity <- partial(
    .f=.plot.specificity,
    "filtering_threshold",
    "Detection threshold (allele frequency)")


plot.num_genes_from_module.sensitivity <- partial(
    .f=.plot.sensitivity,
    "num_genes_from_module",
    "Number of genes sampled per module")


plot.num_genes_from_module.specificity <- partial(
    .f=.plot.specificity,
    "num_genes_from_module",
    "Number of genes sampled per module")


plot.sampling_interval.sensitivity <- partial(
    .f=.plot.sensitivity,
    "sampling_interval",
    "Sampling interval (generations)")


plot.sampling_interval.specificity <- partial(
    .f=.plot.specificity,
    "sampling_interval",
    "Sampling interval (generations)")

plot.timeseries_length.sensitivity <- partial(
    .f=.plot.sensitivity,
    "timeseries_length",
    "Timeseries length (generations)")


plot.timeseries_length.specificity <- partial(
    .f=.plot.specificity,
    "timeseries_length",
    "Timeseries length (generations)")


## This code comes straight from the cowplot documentation:
## https://wilkelab.org/cowplot/articles/plot_grid.html#joint-plot-titles-1
plot.with.title <- function(my.plot, my.title) {

    ## add the title
    title <- ggdraw() + 
        draw_label(
            my.title,
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            ## add margin on the left of the drawing canvas,
            ## so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
        )

    plot.with.title <- plot_grid(
        title, my.plot,
        ncol = 1,
        ## rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
    )
    return(plot.with.title)
}


######################################################################################

## get the power analysis data.
power.analysis.df <- read.csv("../results/SLiM-results/STIMS-power-analysis-data.csv") %>%
    mutate(GeneModule = factor(GeneModule, levels = c("positive", "neutral", "purifying")))

## make confusion matrices.
positive.confusion.matrix <- make.positive.confusion.matrix(power.analysis.df)
negative.confusion.matrix <- make.negative.confusion.matrix(power.analysis.df)

## make positive selection hypermutator plots.
hypermut.positive.filtering_threshold.sensitivity.plot <- positive.confusion.matrix %>%
    plot.filtering_threshold.sensitivity("Hypermutator")

hypermut.positive.filtering_threshold.specificity.plot <- positive.confusion.matrix %>%
    plot.filtering_threshold.specificity("Hypermutator")

hypermut.positive.num_genes_from_module.sensitivity.plot <- positive.confusion.matrix %>%
    plot.num_genes_from_module.sensitivity("Hypermutator")

hypermut.positive.num_genes_from_module.specificity.plot <- positive.confusion.matrix %>%
    plot.num_genes_from_module.specificity("Hypermutator")
    
hypermut.positive.timeseries_length.sensitivity.plot <- positive.confusion.matrix %>%
    plot.timeseries_length.sensitivity("Hypermutator")

hypermut.positive.timeseries_length.specificity.plot <- positive.confusion.matrix %>%
    plot.timeseries_length.specificity("Hypermutator")


hypermut.positive.plots <- plot_grid(
    hypermut.positive.filtering_threshold.sensitivity.plot,
    hypermut.positive.filtering_threshold.specificity.plot,
    hypermut.positive.num_genes_from_module.sensitivity.plot,
    hypermut.positive.num_genes_from_module.specificity.plot,
    hypermut.positive.timeseries_length.sensitivity.plot,
    hypermut.positive.timeseries_length.specificity.plot,
    nrow=3)
                                     

## make negative selection hypermutator plots
hypermut.negative.filtering_threshold.sensitivity.plot <- negative.confusion.matrix %>%
    plot.filtering_threshold.sensitivity("Hypermutator")

hypermut.negative.filtering_threshold.specificity.plot <- negative.confusion.matrix %>%
    plot.filtering_threshold.specificity("Hypermutator")

hypermut.negative.num_genes_from_module.sensitivity.plot <- negative.confusion.matrix %>%
    plot.num_genes_from_module.sensitivity("Hypermutator")

hypermut.negative.num_genes_from_module.specificity.plot <- negative.confusion.matrix %>%
    plot.num_genes_from_module.specificity("Hypermutator")
    
hypermut.negative.timeseries_length.sensitivity.plot <- negative.confusion.matrix %>%
    plot.timeseries_length.sensitivity("Hypermutator")

hypermut.negative.timeseries_length.specificity.plot <- negative.confusion.matrix %>%
    plot.timeseries_length.specificity("Hypermutator")


hypermut.negative.plots <- plot_grid(
    hypermut.negative.filtering_threshold.sensitivity.plot,
    hypermut.negative.filtering_threshold.specificity.plot,
    hypermut.negative.num_genes_from_module.sensitivity.plot,
    hypermut.negative.num_genes_from_module.specificity.plot,
    hypermut.negative.timeseries_length.sensitivity.plot,
    hypermut.negative.timeseries_length.specificity.plot,
    nrow=3)


## make positive selection nonmutator plots.
nonmut.positive.filtering_threshold.sensitivity.plot <- positive.confusion.matrix %>%
    plot.filtering_threshold.sensitivity("Nonmutator")

nonmut.positive.filtering_threshold.specificity.plot <- positive.confusion.matrix %>%
    plot.filtering_threshold.specificity("Nonmutator")

nonmut.positive.num_genes_from_module.sensitivity.plot <- positive.confusion.matrix %>%
    plot.num_genes_from_module.sensitivity("Nonmutator")

nonmut.positive.num_genes_from_module.specificity.plot <- positive.confusion.matrix %>%
    plot.num_genes_from_module.specificity("Nonmutator")
    
nonmut.positive.timeseries_length.sensitivity.plot <- positive.confusion.matrix %>%
    plot.timeseries_length.sensitivity("Nonmutator")

nonmut.positive.timeseries_length.specificity.plot <- positive.confusion.matrix %>%
    plot.timeseries_length.specificity("Nonmutator")


nonmut.positive.plots <- plot_grid(
    nonmut.positive.filtering_threshold.sensitivity.plot,
    nonmut.positive.filtering_threshold.specificity.plot,
    nonmut.positive.num_genes_from_module.sensitivity.plot,
    nonmut.positive.num_genes_from_module.specificity.plot,
    nonmut.positive.timeseries_length.sensitivity.plot,
    nonmut.positive.timeseries_length.specificity.plot,
    nrow=3)


## make negative selection nonmutator plots.
nonmut.negative.filtering_threshold.sensitivity.plot <- negative.confusion.matrix %>%
    plot.filtering_threshold.sensitivity("Nonmutator")

nonmut.negative.filtering_threshold.specificity.plot <- negative.confusion.matrix %>%
    plot.filtering_threshold.specificity("Nonmutator")

nonmut.negative.num_genes_from_module.sensitivity.plot <- negative.confusion.matrix %>%
    plot.num_genes_from_module.sensitivity("Nonmutator")

nonmut.negative.num_genes_from_module.specificity.plot <- negative.confusion.matrix %>%
    plot.num_genes_from_module.specificity("Nonmutator")
    
nonmut.negative.timeseries_length.sensitivity.plot <- negative.confusion.matrix %>%
    plot.timeseries_length.sensitivity("Nonmutator")

nonmut.negative.timeseries_length.specificity.plot <- negative.confusion.matrix %>%
    plot.timeseries_length.specificity("Nonmutator")


nonmut.negative.plots <- plot_grid(
    nonmut.negative.filtering_threshold.sensitivity.plot,
    nonmut.negative.filtering_threshold.specificity.plot,
    nonmut.negative.num_genes_from_module.sensitivity.plot,
    nonmut.negative.num_genes_from_module.specificity.plot,
    nonmut.negative.timeseries_length.sensitivity.plot,
    nonmut.negative.timeseries_length.specificity.plot,
    nrow=3)


hypermut.panelA <- plot.with.title(
    hypermut.positive.plots,
    "STIMS detects positive selection in hypermutator populations")

hypermut.panelB <- plot.with.title(
    hypermut.negative.plots,
    "STIMS detects purifying selection in hypermutator populations")

nonmut.panelA <- plot.with.title(
    nonmut.positive.plots,
    "STIMS detects positive selection in nonmutator populations")

nonmut.panelB <- plot.with.title(
    nonmut.negative.plots,
    "STIMS does not detect purifying selection in nonmutator populations")

hypermut.power.figure <- plot_grid(hypermut.panelA, hypermut.panelB, nrow=2)
ggsave("../results/SLiM-results/hypermutator-power-figure.pdf",
       hypermut.power.figure,
       height = 8)

nonmut.power.figure <- plot_grid(nonmut.panelA, nonmut.panelB, nrow=2)
ggsave("../results/SLiM-results/nonmutator-power-figure.pdf",
       nonmut.power.figure,
       height = 8)
