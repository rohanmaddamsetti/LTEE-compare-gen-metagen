## SLiM-analysis.R by Rohan Maddamsetti.
## This script contains diagnostic code for analyzing SLiM simulation output.
source("metagenomics-library.R")

SLiM.genes <- read.csv("../results/SLiM-results/SLiM_geneIDs.csv")
SLiM.mutation.data <- read.csv("../results/SLiM-results/SLiM-1000gen-FivePercent-Hypermutator.csv")

mut.summary.1 <- SLiM.mutation.data %>%
    group_by(Annotation) %>%
    summarize(count = n())

mut.summary.2 <- SLiM.mutation.data %>%
    group_by(Gene) %>%
    summarize(count = n())

mut.dist.plot <- ggplot(SLiM.mutation.data, aes(x=Gene)) + geom_histogram(stat="count")

non.filtered.mut.data <- read.csv("../results/SLiM-results/SLIM-5000gen-v03.csv")

mut.summary.4 <- old.filtered.SLiM.mutation.data %>%
    group_by(Gene) %>%
    summarize(count = n())

old.filtered.mut.dist.plot <- ggplot(old.filtered.SLiM.mutation.data,
                                     aes(x=Position)) + geom_histogram()
