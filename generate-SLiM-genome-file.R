## generate-SLiM-genome-files.R by Nkrumah Grant and Rohan Maddamsetti
## Create a simple genome.csv file akin to REL606_IDs.csv,
## and create files for each module.

library(tidyverse)

## The genome is divided into 1000 genes, each 1000 bp.
n <- 1000
prefix <- "g"
suffix <- seq(1:n)
Gene <- paste(prefix, suffix, sep="")

end <- c()
start <- c()

for (i in 1:n) {
    end <- c(end, (i*1000))
    start <- c(end - 1000)
    gene_length <- c(end - start)
}

## The first 25 genes have beneficial mutations. 10% of mutations
## in this region are beneficial, and 90% are neutral.
## The next 25 genes have strongly deleterious mutations.
## 40% of mutations in this region are deleterious.
## The rest of the genome is completely neutral.

## neutral, beneficial, deleterious 
## g1 = 0.9, 0.1, 0.0, 
## g2 = 0.6, 0.0, 0.40 
## g3 = 1.0, 0.0, 0.0 

SLiM.genes <- data.frame(Gene, start, end, gene_length) %>%
    ## Let's annotate the genome file.
    mutate(Module = c(rep("Positive", 25),
                      rep("Purifying",25),
                      rep("Neutral", 25),
                      rep("Neutral_Background",925))) %>%
    mutate(PercentageNeutral = c(rep(0.9, 25), rep(0.6, 25),
                                 rep(1.0, 950))) %>%
    mutate(PercentageBeneficial = c(rep(0.1, 25), rep(0.0, 25),
                                    rep(0.0, 950))) %>%
    mutate(PercentageDeleterious = c(rep(0.0, 25), rep(0.4, 25),
                                     rep(0,950)))

write.csv(SLiM.genes, "../results/SLiM-results/SLiM_geneIDs.csv",
          quote = F, row.names = F)

positive.selection.module <- filter(SLiM.genes, Module == "Positive")
purifying.selection.module <- filter(SLiM.genes, Module == "Purifying")
neutral.module <- filter(SLiM.genes, Module == "Neutral")

## write the modules to file.
write.csv(file="../results/SLiM-results/SLiM_positive_module.csv",
          positive.selection.module)
write.csv(file="../results/SLiM-results/SLiM_purifying_module.csv",
          purifying.selection.module)
write.csv(file="../results/SLiM-results/SLiM_neutral_module.csv",
          neutral.module)
