## generate-SLiM-genome-files.R by Nkrumah Grant and Rohan Maddamsetti
## Create a simple genome.csv file akin to REL606_IDs.csv,
## and create files for each module.

library(tidyverse)

## The genome is divided into 4000 genes, each 1000 bp.
n <- 4000
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

## The first 100 genes have beneficial mutations. 10% of mutations
## in this region are beneficial, and 90% are background.
## The next 100 genes have strongly deleterious mutations.
## 40% of mutations in this region are deleterious.
## The next 100 genes are completely neutral.
## The rest of the genome is a weakly deleterious background,
## based on the distribution reported by Lydia Robert et al. (2018).

## beneficial, deleterious, neutral, background
## g1 = 0.1, 0.0, 0.0, 0.9
## g2 = 0.0. 0.4, 0.0, 0.6 
## g3 = 0.0, 0.0, 1.0, 0.0
## g4 = 0.0, 0.0, 0.0, 1.0

SLiM.genes <- data.frame(Gene, start, end, gene_length) %>%
    ## Let's annotate the genome file.
    mutate(Module = c(rep("Positive", 100),
                      rep("Purifying",100),
                      rep("Neutral", 100),
                      rep("Background",3700))) %>%
    mutate(PercentageBeneficial = c(rep(0.1, 100),
                                    rep(0.0, 100),
                                    rep(0.0, 100),
                                    rep(0.0, 3700))) %>%
    mutate(PercentageDeleterious = c(rep(0.0, 100),
                                     rep(0.4, 100),
                                     rep(0.0, 100),
                                     rep(0.0,3700))) %>%
    mutate(PercentageNeutral = c(rep(0.0, 100),
                                 rep(0.0, 100),
                                 rep(1.0, 100),
                                 rep(0.0, 3700))) %>%
    mutate(PercentageBackground = c(rep(0.9, 100),
                                    rep(0.6, 100),
                                    rep(0.0, 100),
                                    rep(1.0, 3700)))

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
