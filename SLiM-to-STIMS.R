## SLiM-to-STIMS.R
## Nkrumah Grant and Rohan Maddamsetti
## take the SliM output file and convert it into input for STIMS.

library(tidyverse)
library(data.table)



## IMPORTANT: we are assuming that each gene is 1000bp long.
annotate.Gene.per.mutation <- function(df, gene.length = 1000) {
    ## annotate the Gene for each mutation.
    df %>%
        mutate(GeneBin = trunc(Position/gene.length)+1) %>% 
        mutate(Gene = paste("g", GeneBin, sep = "")) %>% 
        select(-GeneBin)
}


## This function nicely formats all mutations over time into a dataframe.
## Further processing is needed to generate STIMS input.
SLiM.output.to.dataframe <- function(SLiM.outfile, pop.name,
                                     freq_threshold=0.01, Ne=1e5) {
    ## all mutations in the population, sampled every 100 generations.
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

## This function makes STIMS input-- only one mutation per row,
## representing the first time the mutation was observed in the population.
SLiM.output.to.STIMS.input <- function(SLiM.output, pop.name,
                                       freq_threshold=0.01, Ne=1e5) {
    SLiM.output.to.dataframe(SLiM.output, pop.name, freq_threshold, Ne) %>%
        ## IMPORTANT: There can only be one row per mutation.
        group_by(ID, Annotation, Position, Gene, Population, t0) %>%
        arrange(Generation) %>%
        ## only take the first row (the first generation at which allele_freq > threshold)
        filter(row_number() == 1) %>%
        arrange(ID) ## to double-check that each ID is unique    
}

## Simulated Hypermutator data.
d1 <- SLiM.output.to.STIMS.input(
    "../results/SLiM-results/SLiM-replicate-runs/SLiM_Ne1000000_mu10-8_numgens5000_Rep999.txt",
##    "../results/SLiM-results/SLiM_Ne1000000_mu10-8_numgens5000.txt",
    "Hypermutator",
    freq_threshold = 0.01, Ne = 1e6)
## write out the STIMS input file.
write.csv(
    d1,
    "../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator_Rep999.csv",
##    "../results/SLiM-results/SLiM-5000gen-OnePercent-Hypermutator.csv",
    quote = F, row.names = F)

## Simulated Nonmutator data.
##d2 <- SLiM.output.to.STIMS.input(
##    "../results/SLiM-results/SLiM_Ne1000000_mu10-10_numgens5000.txt",
##    "Nonmutator",
##    freq_threshold = 0.01, Ne = 1e6)
#### write out the STIMS input file.
##write.csv(
##    d2,
##    "../results/SLiM-results/SLiM-5000gen-OnePercent-Nonmutator.csv",
##    quote = F, row.names = F)

