## analyze-STIMS-power.R by Rohan Maddamsetti and Nkrumah Grant. 
## This script takes the output of generate-STIMS-power-data.jl,
## and makes figures.

library(tidyverse)
library(cowplot)
library(stats)

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


#Here, I calculate the false discovery rate using the list of, 
#calculated p-values from the simulation data. I am not sure what, 
#to do with the neutral p-values.

#TODO 
#Make sure to calculate Q.val for neutral case. As of now, the values 
#in the table are the p-values. I think a range between the upper and, 
#lower bound should be reported for the neutral case. Discuss with
#Rohan. 

calculate.Qval <- function(df) {
  df %>% 
    group_by(MutatorStatus, GeneModule, Variable, VariableValue) %>%  
    summarise(Qval = p.adjust(Pval, method = "fdr", n = length(Pval))) %>% 
    group_by(MutatorStatus, GeneModule, Variable, VariableValue) %>% 
    summarise(avg.Qval = mean(Qval), n = n()) %>% 
    mutate(avg.false.positives = ifelse(GeneModule == "purifying", abs(1-avg.Qval), avg.Qval)) %>% 
    mutate(avg.false.positives = ifelse(GeneModule == "positive", abs(0-avg.Qval), avg.false.positives))
}

df1 <- calculate.Qval(timeseries.length.df)
df2 <- calculate.Qval(sampling.interval.df)
df3 <- calculate.Qval(filtering.threshold.df)
df4 <- calculate.Qval(num.genes.from.module.df)

write.csv(false.discovery.rate <- rbind(df1, df2, df3, df4), "../results/SliM-results/false-discovery-rate.csv", row.names = F)






