## analyze-STIMS-power.R by Rohan Maddamsetti.
## This script takes the output of generate-STIMS-power-data.jl,
## and makes figures.

library(tidyverse)
library(cowplot)

power.analysis.df <- read.csv("../results/SLiM-results/STIMS-power-analysis-data.csv") %>%
    mutate(GeneModule = factor(GeneModule, levels = c("positive", "neutral", "purifying")))

## now, let's make the figures!
timeseries.length.df <- power.analysis.df %>%
    filter(Variable == "timeseries_length")

sampling.interval.df <- power.analysis.df %>%
    filter(Variable == "sampling_interval")

filtering.threshold.df <- power.analysis.df %>%
    filter(Variable == "filtering_threshold")

num.genes.from.module.df <- power.analysis.df %>%
    filter(Variable == "num_genes_from_module")

make.power.plot <- function(df) {
    df %>%
        ggplot(aes(x = VariableValue, y = Pval)) +
        theme_classic() +
        geom_point() +
        facet_grid(MutatorStatus~GeneModule) +
        ylab("p-value") +
        geom_hline(yintercept = 0.05,linetype="dashed", color="red") +
        geom_hline(yintercept = 0.95,linetype="dashed", color="red")
}

timeseries.plot <- make.power.plot(timeseries.length.df) +
    ggtitle("Length of Simulated Evolution Experiment") +
    xlab("Generations")

sampling.plot <- make.power.plot(sampling.interval.df) +
    ggtitle("Interval between Samples") +
    xlab("Generations")

filtering.threshold.plot <- make.power.plot(filtering.threshold.df) +
    ggtitle("Metagenomic Allele Frequency Detection Limit") +
    xlab("Allele Frequency Threshold")

num.genes.from.module.plot <- make.power.plot(num.genes.from.module.df) +
    ggtitle("Number of genes sampled from ground truth module") +
    xlab("Number of genes per sample")


power.figure <- plot_grid(timeseries.plot,
                          sampling.plot,
                          filtering.threshold.plot,
                          num.genes.from.module.plot,
                          labels = c("A","B","C","D"))
ggsave("../results/SLiM-results/STIMS-power.pdf")
