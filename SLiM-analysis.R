## SLiM-analysis.R by Rohan Maddamsetti.
## This script contains diagnostic code for analyzing SLiM simulation output.

library(tidyverse)

SLiM.genes <- read.csv("../results/SLiM-results/SLiM_geneIDs.csv")
SLiM.mutation.data <- read.csv(
    "../results/SLiM-results/SLiM-1000gen-FivePercent-Hypermutator.csv") %>%
    ## add gene metadata.
    left_join(SLiM.genes)

mut.summary.1 <- SLiM.mutation.data %>%
    group_by(Annotation) %>%
    summarize(count = n())

mut.summary.2 <- SLiM.mutation.data %>%
    group_by(Gene) %>%
    summarize(count = n())

## Let's plot the distribution of mutations over the genome.
mut.dist.plot <- ggplot(SLiM.mutation.data, aes(x=Gene)) + geom_histogram(stat="count")
ggsave("../results/SLiM-results/SLiM-genomic-mutation-distribution.pdf",
       mut.dist.plot)


calc.SLiM.gene.mutation.density <- function(SLiM.mutation.data, SLiM.genes) {
    ## Examine the distribution of various classes of mutations across genes in the
    ## SLiM metagenomics data. Which genes are enriched? Which genes are depleted?
    ## IMPORTANT: GENES WITH NO MUTATIONS ARE OMITTED.
    density.df <- SLiM.mutation.data %>%
        mutate(Gene=as.factor(Gene)) %>%
        group_by(Gene,gene_length) %>%
        summarize(mut.count=n()) %>%
        ungroup() %>%
        mutate(mut.density=mut.count/gene_length) %>%
        arrange(desc(mut.density))
    return(density.df)
}


make.SLiM.gene.all.mut.density.df <- function(gene.mutation.data, SLiM.genes) {
    gene.mutation.density <- calc.SLiM.gene.mutation.density(gene.mutation.data)

    ## This is to include zeros.
    gene.mutation.densities <- SLiM.genes %>%
        full_join(gene.mutation.density)

    ## CRITICAL STEP: replace NAs with zeros.
    ## We need to keep track of genes that haven't been hit by any mutations
    gene.mutation.densities[is.na(gene.mutation.densities)] <- 0
    gene.mutation.densities <- as_tibble(gene.mutation.densities)

    return(gene.mutation.densities)
}

gene.mut.density.df <- make.SLiM.gene.all.mut.density.df(
    SLiM.mutation.data, SLiM.genes) %>%
    arrange(desc(mut.density)) %>%
    mutate(Gene_ranked_by_mut_density = row_number())

## Let's plot the ranked density of mutations per gene.
mut.density.plot <- gene.mut.density.df %>%
    ggplot(aes(x=Gene_ranked_by_mut_density, y=mut.density)) +
    theme_classic() + geom_point()

write.csv(gene.mut.density.df,
          "../results/SLiM-results/SLiM-gene-mutation-density.csv")
ggsave("../results/SLiM-results/SLiM-gene-mutation-density.pdf",
       mut.density.plot)
