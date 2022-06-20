## SLiM-plot-mutation-trajectories.R
## Nkrumah Grant and Rohan Maddamsetti
## Take SLiM output data and create a mutation trajectory plot.

library(tidyverse)
library(cowplot)
library(gghighlight)


## IMPORTANT: we are assuming that each gene is 1000bp long.
annotate.Gene.per.mutation <- function(df, gene.length = 1000) {
    ## annotate the Gene for each mutation.
    df %>%
        mutate(GeneBin = trunc(Position/gene.length)+1) %>% 
        mutate(Gene = paste("g", GeneBin, sep = "")) %>% 
        select(-GeneBin)
}


SLiM.output.to.dataframe <- function(SLiM.outfile, pop.name,
                                     freq_threshold=0.01, Ne=1e5) {
    SLiM.outfile %>%
        data.table::fread(header = F, sep = " ") %>%
        ## Remove unnecessary columns from SLiM output. 
        select(c("V2", "V5" , "V6", "V7", "V10", "V11", "V12")) %>%
        rename(Generation = V2) %>%
        rename(ID = V5) %>%
        rename(Annotation = V6) %>%
        rename(Position = V7) %>%
        rename(Population = V10) %>%
        rename(t0 = V11) %>%
        rename(prevalence = V12) %>%
        mutate(Annotation = recode(Annotation,
                                   m1 = "beneficial",
                                   m2 = "deleterious",
                                   m3 = "neutral",
                                   m4 = "background")) %>%
        mutate(Population = recode(Population, p1 = pop.name)) %>%
        mutate(Position = as.numeric(Position)) %>%
        ## annotate the Gene for each mutation.
        annotate.Gene.per.mutation() %>%
        ## Convert prevalence to allele frequency.
        ## The denominator is the number of individuals in the SliMulation. 
        mutate("allele_freq" = prevalence/Ne) %>%
        ## Filter mutations with allele frequencies above the sampling threshold.
        filter(allele_freq > freq_threshold) %>%
        arrange(ID, Generation) ## arrange each mutation by generation.
}


#Make mutation trajectory plots
make.mutation.trajectory.plot <- function(df) {
    ggplot(df, aes(x=Generation, y=allele_freq, colour = as.factor(ID)))+
        geom_line() +
        gghighlight(max(allele_freq) == 1) +
        xlab("Generations") +
        ylab("Allele frequency") +
        theme_classic() +
        theme(axis.text.x = element_text(colour = "black", size = 14,
                                         margin = (margin(t = 5, b=5)))) +
        theme(axis.text.y = element_text(colour = "black", size = 14,
                                         margin = (margin(l = 5, r=5)))) +
        theme(axis.title.x = element_text(size = 14)) +
        theme(axis.title.y = element_text(size = 14)) +
        theme(panel.border = element_blank(), axis.line = element_line()) +
        theme(legend.position = "none")
}


## all mutations in the population, sampled every 100 generations.
## IMPORTANT: Pass in the correct Ne and filtering threshold.
hypermut.df <- SLiM.output.to.dataframe(
    "../results/SLiM-results/SLiM_Ne1000000_mu10-8_numgens5000.txt",
    "Hypermutator",
    freq_threshold = 0.00001, ## 1e-5 precision for nice curves.
    Ne = 1e6)

nonmut.df <- SLiM.output.to.dataframe(
    "../results/SLiM-results/SLiM_Ne1000000_mu10-10_numgens5000.txt",
    "Nonmutator",
    freq_threshold = 0.0000001, ## 1e-7 precision for nice curves.
    Ne = 1e6)

hypermut.p <- make.mutation.trajectory.plot(hypermut.df)
ggsave("../results/SLiM-results/figures/hypermut-trajectory.pdf", hypermut.p, height=2.5, width=6.5)

nonmut.p <- make.mutation.trajectory.plot(nonmut.df) 
ggsave("../results/SLiM-results/figures/nonmut-trajectory.pdf", nonmut.p, height=2.5, width=6.5)
