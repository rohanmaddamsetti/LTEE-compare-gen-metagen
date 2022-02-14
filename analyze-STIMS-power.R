## analyze-STIMS-power.R by Rohan Maddamsetti.
## This script takes the output of generate-STIMS-power-data.jl,
## and makes figures.

library(tidyverse)

power.analysis.df <- read.csv("../results/SLiM-results/STIMS-power-analysis-data.csv")

## now, let's make the figures!
##timeseries_length_df = rsubset(power_analysis)

##time_series_plot = ggplot(timeseries_length_df, aes(x = :VariableValue, y = :Pval)) +
##    theme_classic() +
##    geom_boxplot() +
##    facet_grid(R"MutatorStatus~GeneModule") +
##    ggtitle("Length of Simulated Evolution Experiment") +
##    xlab("Generations")
