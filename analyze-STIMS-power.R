## analyze-STIMS-power.R by Rohan Maddamsetti and Nkrumah Grant. 
## This script takes the output of generate-STIMS-power-data.jl,
## and makes figures.

library(tidyverse)
library(cowplot)
library(boot)


calculate.Qval <- function(df) {
    ## Calculate FDR-corrected p-values for each treatment, and average.
    df %>%
        group_by(MutatorStatus, GeneModule, Variable, VariableValue) %>%  
        summarise(Qval = stats::p.adjust(Pval, method = "fdr", n = n())) %>% 
        group_by(MutatorStatus, GeneModule, Variable, VariableValue) %>% 
        summarise(avg.Qval = mean(Qval), n = n()) %>%
        ## take the complement (1-q) for purifying selection.
        mutate(avg.Qval = ifelse(GeneModule == "purifying", 1 - avg.Qval, avg.Qval))
}


calc.bootstrap.conf.int <- function(vec) {
  ## bootstrap confidence intervals around the mean.
  ## Use the boot package to calculate fancy BCA intervals.
  
  mean.boot <- function(x,ind) {
      ## define this type of mean-calculating function to pass to the boot function.
      ## from: https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/9-18.pdf
      ## as it seems this link has been taken down, I downloaded these notes
      ## using the Wayback Machine from the Internet Archive.
      return(c(mean(x[ind]), var(x[ind])/length(x)))
  }

    Nbootstraps <- 10000

    ## remove any NA values from vec.
    vec <- vec[!is.na(vec)]
    
    if (length(vec)) {    
        out <- boot(vec,mean.boot,Nbootstraps)
        ## handle bad inputs in bootstrapping confidence intervals
        ci.result <- tryCatch(boot.ci(out,type="bca"), error= function(c) return(NULL))
        if (is.null(ci.result)) {
            Left <- NA
            Right <- NA
        } else {
            Left <- ci.result$bca[4]
            Right <- ci.result$bca[5]
        }
    } else { ## vec only contained NA values.
        Left <- NA
        Right <- NA
    }   
    final <- c(Left, Right)
    return(final)
}


calc.Pval.confint.df <- function(df.group) {
    Pval.conf.int <- calc.bootstrap.conf.int(df.group$Pval)
    results <- data.frame(MutatorStatus = unique(df.group$MutatorStatus),
                          GeneModule = unique(df.group$GeneModule),
                          Variable = unique(df.group$Variable),
                          VariableValue = unique(df.group$VariableValue),
                          mean_Pval = mean(df.group$Pval),
                          Left = Pval.conf.int[1],
                          Right = Pval.conf.int[2])
    return(results)
}


make.pval.panel <- function(panel.df) {
    panel.df %>%
        ## take the complement for purifying selection p-values.
        mutate(Pval = ifelse(GeneModule == "purifying", 1 - Pval, Pval)) %>%
        group_split(MutatorStatus, GeneModule, Variable, VariableValue) %>%
        map_dfr(.f=calc.Pval.confint.df) %>%
        ## for these particular data, we get NAs
        ## when the p-value is 0, with no variance.
        replace_na(list(Left = 0,Right=0)) %>% 
        ggplot(aes(x = VariableValue, y = mean_Pval)) +
        theme_classic() +
        geom_point() +
        ## plot CIs.
        geom_errorbar(aes(ymin=Left,ymax=Right)) +
        facet_grid(MutatorStatus~GeneModule) +
        ylab("p-value") +
        geom_hline(yintercept = 0.05,linetype="dashed", color="red")
}


make.pval.plot <- function(plot.df) {

    timeseries.length.df <- filter(plot.df, Variable == "timeseries_length")
    sampling.interval.df <- filter(plot.df, Variable == "sampling_interval")
    filtering.threshold.df <- filter(plot.df, Variable == "filtering_threshold")
    num.genes.from.module.df <- filter(plot.df, Variable == "num_genes_from_module")

    timeseries.length.panel <- make.pval.panel(timeseries.length.df) +
        ggtitle("Length of Evolution Experiment") +
        xlab("Generations")

    sampling.interval.panel <- make.pval.panel(sampling.interval.df) +
        ggtitle("Interval between Samples") +
        xlab("Generations")

    filtering.threshold.panel <- make.pval.panel(filtering.threshold.df) +
        ggtitle("Allele Frequency Detection Limit") +
        xlab("Allele Frequency Threshold")

    num.genes.from.module.panel <- make.pval.panel(num.genes.from.module.df) +
        ggtitle("Number of genes per random module") +
        xlab("Number of genes per sample")
    

    pval.plot <- plot_grid(timeseries.length.panel,
                           sampling.interval.panel,
                           filtering.threshold.panel,
                           num.genes.from.module.panel)
    return(pval.plot)
}

###########################################################################################
## get the power analysis data.
power.analysis.df <- read.csv("../results/SLiM-results/STIMS-power-analysis-data.csv") %>%
    mutate(GeneModule = factor(GeneModule, levels = c("positive", "neutral", "purifying")))

###########################################################################################
## for the p-value plots, show Hypermutators and Nonmutators separately,

positive.hypermut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Hypermutator") %>%
    filter(GeneModule == "positive") %>%
    make.pval.plot()

neutral.hypermut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Hypermutator") %>%
    filter(GeneModule == "neutral") %>%
    make.pval.plot()

purifying.hypermut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Hypermutator") %>%
    filter(GeneModule == "purifying") %>%
    make.pval.plot()

positive.nonmut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Nonmutator") %>%
    filter(GeneModule == "positive") %>%
    make.pval.plot()

neutral.nonmut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Nonmutator") %>%
    filter(GeneModule == "neutral") %>%
    make.pval.plot()

    purifying.nonmut.plot <- power.analysis.df %>%
    filter(MutatorStatus == "Nonmutator") %>%
    filter(GeneModule == "purifying") %>%
    make.pval.plot()

ggsave("../results/SLiM-results/positive-hypermut-plot.pdf", positive.hypermut.plot)
ggsave("../results/SLiM-results/neutral-hypermut-plot.pdf", neutral.hypermut.plot)
ggsave("../results/SLiM-results/purifying-hypermut-plot.pdf", purifying.hypermut.plot)

ggsave("../results/SLiM-results/positive-nonmut-plot.pdf", positive.nonmut.plot)
ggsave("../results/SLiM-results/neutral-nonmut-plot.pdf", neutral.nonmut.plot)
ggsave("../results/SLiM-results/purifying-nonmut-plot.pdf", purifying.nonmut.plot)
####################################################################





