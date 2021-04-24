## Poisson-selection-analysis.R by 
## Nkrumah Grant
## https://docs.rstudio.com/tutorials/user/using-python-with-rstudio-and-reticulate/

## TODO: fix the bug that Nkrumah found in upstream annotation code:
## Noticed the gene length was listed as "0" for ydhP. Actual gene length is 1167.
## I added a step to manually add the value. 

library(tidyverse)

#Import all of the data files that are necessary for the analysis

REL606.genes <- read.csv("../results/REL606_IDs.csv", as.is = T)
meta.genomes <- read.csv("../results/LTEE-metagenome-mutations.csv", as.is =T)
mod.regulators <- read.csv("../data/rohans-I-modulons-to-regulators.csv", as.is =T)
nmut.parallelism <- read.csv("../data/tenaillon2016-nonmutator-parallelism.csv", as.is=T)
mut.parallelism <- read.csv("../data/tenaillon2016-mutator-parallelism.csv", as.is=T)
couce <- read.csv("../data/Couce2017-LTEE-essential.csv")

## count number of mutated genes in each of the populations
## and then calculate the overall density of mutations
## Length of REL606 Genome = 4,629,812 bp (Jeong et al. 2009)
Sum.mutations <- meta.genomes %>%
  group_by(Population) %>%
  summarise(sum.mutations = n(), background.mut.density = sum.mutations/4629812)  

## Count the total number of mutations that occur in each gene.
## I can filter the synonymous mutations out of this analysis if necessary.
Pois.data.analysis <- meta.genomes %>% group_by(Population) %>%
    count(Gene) %>% rename(mut.count = n) %>%
    left_join(select(mut.parallelism, Gene.name, Coding.length), by = c("Gene" = "Gene.name")) 

## Noticed the gene length was listed as "0" for ydhP. Actual gene length is 1167.
## I added a step to manually add the value. 
## manually add gene length for ydhP
Pois.data.analysis$Coding.length[Pois.data.analysis$Coding.length == '0'] <- as.numeric('1167')

Pois.data.analysis <- Pois.data.analysis %>%
  left_join(select(Sum.mutations, Population, background.mut.density)) %>% 
    mutate(pois.lambda = background.mut.density * Coding.length)

## calculate poisson probability:
## Here I hard-coded the poisson distribution as presented in the Kinnersley paper.
## One can also use the dpois function in base R. The results are the same. 
Pois.data.analysis <- Pois.data.analysis %>%
    mutate(my.pois.pval = (pois.lambda^mut.count * exp(-pois.lambda)/factorial(mut.count))) %>%
    mutate(pois.pval = dpois(x = mut.count,lambda = pois.lambda))

## assert that the numbers are almost equal per numerical precision.
stopifnot(all.equal(Pois.data.analysis$my.pois.pval, Pois.data.analysis$pois.pval))

## write results to file.
write.csv(Pois.data.analysis, file = "../results/gene-modules/poisson.data.analysis.csv")


#################################

## Rohan messing around with results here.

significant.genes <- filter(Pois.data.analysis, pois.pval < 0.001)

## how many of these are in the nonmutator parallelism table?
## we have to filter on the top G-scoring genes first.

top.nmut.parallelism <- nmut.parallelism %>% filter(Observed.nonsynonymous.mutation > 1)

nmut.parallel.significant.genes <- significant.genes %>%
    filter(Gene %in% top.nmut.parallelism$Gene.name)
length(unique(nmut.parallel.significant.genes$Gene)) ## 22 genes total.

## what about in Supplementary Table 3 of Good et al. (2017)?
good.S3.table <- read.csv("../data/Good2017-TableS3.csv")

good.S3.significant.genes <- significant.genes %>%
    filter(Gene %in% good.S3.table$Gene)
length(unique(good.S3.significant.genes$Gene)) ## 48 genes total.

## Interesting! results from non-mutators genomes and Good et al. analysis
## do not overlap as much as I expected.
top.nmut.parallelism$Gene.name %in% good.S3.table$Gene

## let's look at the genes which are in neither of the previous reported results.

new.significant.genes <- significant.genes %>%
    filter(!(Gene %in% top.nmut.parallelism$Gene.name)) %>%
    filter(!(Gene %in% good.S3.significant.genes$Gene))

## cpsG, sulA, yeiB, yieN, ykgE, yobF
parallel.new.significant.genes <- new.significant.genes %>% group_by(Gene) %>%
    summarize(count = n()) %>% filter(count > 1)

new.significant.genes %>%
    filter(Gene %in% parallel.new.significant.genes$Gene) %>%
    data.frame()

meta.genomes %>% filter(Population == "Ara-4") %>%
    filter(Gene == "flu")

## pretty cool: quite a bit of parallel evolution here.
parallel.bp.muts <- meta.genomes %>% filter(Gene %in% new.significant.genes$Gene) %>%
    group_by(Position, Gene) %>%
    summarize(count = n()) %>% filter(count > 1)

## more super cool bp-level parallel evolution.
parallel.bp.muts2 <- meta.genomes %>% filter(Gene %in% significant.genes$Gene) %>%
    group_by(Position, Gene) %>%
    summarize(count = n()) %>% filter(count > 1)

## wild-- even more bp-level parallel evolution in these new parallel evolution genes.
parallel.bp.muts3 <- meta.genomes %>%
    filter(Gene %in% parallel.new.significant.genes$Gene) %>%
    group_by(Position, Gene) %>%
    summarize(count = n()) %>% filter(count > 1)

## one thing to keep in mind: some of this stuff is real, some parallelism seems
## to be phase-variation at hypermutable contingency loci. Still could be selection!
## will probably have to work through, case by case.
## may also be some cases caused by things like recombination/gene conversion as well,
## or phenomena that "look" like enrichment of mutations due to parallel evolution or
## selection, but caused by some kind of more complicated evolutionary scenario.

## POTENTIAL TODO: more sophisticated/complicated background mutation density for
## poisson model. Could fit a local mutation density by binning the genome, like
## in divergent mutation bias paper, and some of the STIMS code allowing for
## regional mutation variation, etc. May not be worth the extra complexity though.
