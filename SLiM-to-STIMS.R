## SLiM-to-STIMS.R
## Nkrumah Grant and Rohan Maddamsetti
## take the SliM output file and convert it into input for STIMS.

library(tidyverse)
library(gghighlight)
library(cowplot)
library(magick)
library(ggplotify)

freq_threshold <- 0.05

## This file reflects all mutations that were segregating
## in the population as were sampled every 100 
## generations.
d <- read.csv("../results/SLiM-results/20211212-Mutator.txt", header = F, sep = " ") %>%
    ## Remove unnecessary columns from SLiM output. 
    select(c("V2", "V5" , "V6", "V7", "V10", "V11", "V12")) %>%
    rename(Generation = V2) %>%
    rename(ID = V5) %>%
    rename(Annotation = V6) %>%
    rename(Position = V7) %>%
    rename(Population = V10) %>%
    rename(t0 = V11) %>%
    rename(prevalence = V12) %>%
    mutate(Annotation = recode(Annotation, m1 = "neutral",
                               m2 = "beneficial", m3 = "deleterious")) %>%
    mutate(Population = recode(Population, p1 = "Hypermutator")) %>%
    mutate(Position = as.numeric(Position)) %>%
    ## annotate the Gene for each mutation.
    mutate(P1000 = trunc(Position/1000)+1) %>% 
    mutate(Gene = paste("g", P1000, sep = "")) %>% 
    select(-P1000) %>%
    ## Create new column that converts prevalence to allele frequency.
    ## The denominator is the number of individuals in the SliMulation. 
    mutate("allele_freq" = prevalence/5e5) %>%
    ## Filter mutations with allele frequencies above 5% threshold. 
    ## Maybe look at 5% and 10% frequency thresholds. 
    filter(allele_freq > freq_threshold) %>%
    ## IMPORTANT: There can only be one row per mutation.
    group_by(ID, Annotation, Position, Population, t0) %>%
    arrange(Generation) %>%
    ## only take the first row (the first generation at which allele_freq > threshold)
    filter(row_number() == 1) %>%
    arrange(ID) ## to double-check that each ID is unique

## write out the STIMS input file.
write.csv(d, "../results/SLiM-results/SLiM-1000gen-FivePercent-Hypermutator.csv",
          quote = F, row.names = F)
