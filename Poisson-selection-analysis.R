## Poisson-selection-analysis.R by 
## Nkrumah Grant and Rohan Maddamsetti.

source("metagenomics-library.R")

## get the lengths of all protein-coding genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
## Do by running:
## python printEcoliIDs.py -i ../data/REL606.7.gbk --mask ../data/REL606.L20.G15.P0.M35.RM-edited.mask.gd > ../results/REL606_IDs.csv
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
    mutate(gene_length=strtoi(gene_length)) %>%
    mutate(oriC_start=rotate.REL606.chr(start,"oriC")) %>%
    mutate(oriC_end=rotate.REL606.chr(end,"oriC"))

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")

LTEE_POPS <- length(nonmutator.pops) + length(hypermutator.pops)

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
gene.mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    ## important: inner_join only includes genes that pass REL606 annotation
    ## filters. In particular, genes that overlap with masked regions are
    ## excluded, even if the metagenomic data included mutations that
    ## are called in non-repetitive regions of those genes.
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

#########################################################################
## Kinnersley et al. (2021) style analysis, using background density.

## count number of mutated genes in each of the populations
## and then calculate the overall density of mutations
## since some regions are masked in the metagenomic data,
## use the total length of genes in REL606.genes as the normalizing factor.
TOTAL_GENE_BP <- sum(REL606.genes$gene_length)

pop.sum.mutations <- gene.mutation.data %>%
  group_by(Population) %>%
    summarise(pop.mutation.total = n(),
              background.mut.density = pop.mutation.total/TOTAL_GENE_BP)

pop.poisson.data.analysis <- gene.mutation.data %>%
    group_by(Population,Gene,locus_tag,gene_length) %>%
    summarize(mut.count = n()) %>% ## keep empty groups!
    ungroup() %>%
    ## to include genes with zero mutations in the analysis, join to the
    ## cartesian product of REL606.genes and all twelve populations.
    right_join(crossing(REL606.genes, ## cartesian product here
                        distinct(select(gene.mutation.data,Population)))) %>%
    mutate(mut.count = replace_na(mut.count, 0)) %>% ## put in the zeros.
    left_join(pop.sum.mutations) %>% ## add in the background mut density
    mutate(pois.lambda = background.mut.density * gene_length) %>%
    ## calculate poisson p-value:
    mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
    mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                   lower.tail=FALSE)) %>%
    ## control for the number of statistical comparisons being made:
    mutate(fdr.corrected.lower.pois.prob = p.adjust(
               lower.pois.prob, method = "fdr")) %>%
    mutate(fdr.corrected.upper.pois.prob = p.adjust(
               upper.pois.prob, method = "fdr")) %>%
    mutate(min.tail.fdr.corrected.pois.prob = pmin(
               fdr.corrected.lower.pois.prob,
               fdr.corrected.upper.pois.prob)) %>%
    filter(min.tail.fdr.corrected.pois.prob < 0.01) %>%
    select(Population, Gene, locus_tag, gene_length, mut.count, product,
           fdr.corrected.lower.pois.prob,
           fdr.corrected.upper.pois.prob,
           min.tail.fdr.corrected.pois.prob) %>%
    arrange(Population, min.tail.fdr.corrected.pois.prob)

write.csv(pop.poisson.data.analysis, "../results/FileS3_B.csv")

################################################################################
## summary over all LTEE.
gene.mutation.data.total <- gene.mutation.data %>%
        mutate(mutation.total = nrow(.)) %>%
        mutate(background.mut.density = mutation.total/TOTAL_GENE_BP) %>%
        select(mutation.total,background.mut.density) %>% distinct()

all.ltee.poisson.data.analysis <- gene.mutation.data %>%
    group_by(Gene,locus_tag,gene_length) %>%
    ## filter out synonymous mutations.
    filter(Annotation != "synonymous") %>%
    summarize(mut.count = n()) %>% ## keep empty groups!
    ungroup() %>%
    ## include genes with zero mutations by joining to REL606.genes.
    right_join(REL606.genes) %>%
    mutate(mut.count = replace_na(mut.count, 0)) %>% ## put in the zeros.
    ## add in the background mut density
    left_join(gene.mutation.data.total, by = character()) %>% 
    mutate(pois.lambda = background.mut.density * gene_length) %>%
    ## calculate poisson p-value:
    mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
    mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                   lower.tail=FALSE)) %>%
    ## control for the number of statistical comparisons being made:
    mutate(fdr.corrected.lower.pois.prob = p.adjust(
               lower.pois.prob, method = "fdr")) %>%
    mutate(fdr.corrected.upper.pois.prob = p.adjust(
               upper.pois.prob, method = "fdr")) %>%
    mutate(min.tail.fdr.corrected.pois.prob = pmin(
               fdr.corrected.lower.pois.prob,
               fdr.corrected.upper.pois.prob)) %>%
    filter(min.tail.fdr.corrected.pois.prob < 0.01) %>%
    select(Gene, locus_tag, gene_length, mut.count, product,
           fdr.corrected.lower.pois.prob,
           fdr.corrected.upper.pois.prob,
           min.tail.fdr.corrected.pois.prob) %>%
    arrange(fdr.corrected.upper.pois.prob)

write.csv(all.ltee.poisson.data.analysis, "../results/FileS3_A.csv")

################################################################################

All.LTEE.Poisson.module.analysis <- function(gene.vec, gene.mutation.data,
                                             REL606.genes) {
    ## do the Poisson analysis of modules in STIMS here.
    ## all mutations, including synonymous mutations, are included.

    TOTAL_GENE_BP <- sum(REL606.genes$gene_length)
    
    gene.mutation.data.total <- gene.mutation.data %>%
        mutate(mutation.total = nrow(.)) %>%
        mutate(background.mut.density = mutation.total/TOTAL_GENE_BP) %>%
        select(mutation.total,background.mut.density) %>% distinct()
    
    gene_vec_length <- sum(filter(REL606.genes,Gene %in% gene.vec)$gene_length)
    
    gene.mutation.data %>%
        ## filter for genes in gene.vec.
        filter(Gene %in% gene.vec) %>%
        mutate(group_col="all") %>%
        ## group all data pertaining to the module of interest.
        group_by(group_col) %>%
        summarize(mut.count = n()) %>% ## keep empty groups!
        ungroup() %>%
        mutate(mut.count = replace_na(mut.count, 0)) %>% ## put in the zeros.
        ## add in the background mut density.
        left_join(gene.mutation.data.total, by = character()) %>%
        ## CRITICAL LINE: normalize by the length of all genes in gen.vec.
        mutate(pois.lambda = background.mut.density * gene_vec_length) %>%
        ## calculate poisson p-value:
        mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
        mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                       lower.tail=FALSE)) %>%
        ## control for the number of statistical comparisons being made:
        mutate(fdr.corrected.lower.pois.prob = p.adjust(
                   lower.pois.prob, method = "fdr")) %>%
        mutate(fdr.corrected.upper.pois.prob = p.adjust(
                   upper.pois.prob, method = "fdr")) %>%
        mutate(min.tail.fdr.corrected.pois.prob = pmin(
                   fdr.corrected.lower.pois.prob,
                   fdr.corrected.upper.pois.prob)) %>%
        select(background.mut.density,pois.lambda,fdr.corrected.lower.pois.prob,
               fdr.corrected.upper.pois.prob,min.tail.fdr.corrected.pois.prob) %>%
        arrange(min.tail.fdr.corrected.pois.prob)
}

## This is the function that I will use throughout,
##since the last two arguments are always the same.
Run.All.LTEE.Poisson.module.analysis <- partial(
    All.LTEE.Poisson.module.analysis,
    gene.mutation.data=gene.mutation.data,
    REL606.genes=REL606.genes, ... = )

neutral.genes <- read.csv("../data/neutral_compilation.csv",
                          header=TRUE, as.is=TRUE) %>%
    ## only analyze loci that pass filters.
    filter(Gene %in% REL606.genes$Gene)

neutral.gene.result <- Run.All.LTEE.Poisson.module.analysis(neutral.genes$Gene) %>%
    select(-background.mut.density, -pois.lambda)

## Note: this is a dataframe with one column, called regulator.
Imodulon.regulators <- read.csv(
    "../data/rohans-I-modulons-to-regulators.csv", as.is =T) %>%
    select(regulator) %>% drop_na() %>% distinct()
    
Imodulon.regulator.result <- Run.All.LTEE.Poisson.module.analysis(
    Imodulon.regulators$regulator) %>%
    select(-background.mut.density, -pois.lambda)

couce.essential.genes <- read.csv("../data/Couce2017-LTEE-essential.csv")

couce.essential.result <- Run.All.LTEE.Poisson.module.analysis(
    couce.essential.genes$Gene) %>%
    select(-background.mut.density, -pois.lambda)

## get mutation parallelism in the LTEE genomes published in
## Tenaillon et al. (2016).
## This is used for filtering essential genes
## in the purifying selection control analysis,
## and is used in the positive selection control analysis.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv') %>% ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)

hypermut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv') %>%
    ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)    

## top genes in non-mutators.
top.nonmut.genomics <- top_n(nonmut.genomics, 50, wt=G.score)
## top genes in hypermutators.
top.hypermut.genomics <- top_n(hypermut.genomics, 50, wt=G.score)

top.nonmut.result <- Run.All.LTEE.Poisson.module.analysis(
    top.nonmut.genomics$Gene.name)

top.hypermut.result <- Run.All.LTEE.Poisson.module.analysis(
    top.hypermut.genomics$Gene.name)

## core genes
core.gene.assignments <- read.csv("../data/Maddamsetti2017-core-summary.csv")
core.genes <- core.gene.assignments %>% filter(panortholog == TRUE) %>%
    inner_join(REL606.genes) %>% ## make sure these genes passed the filters
    filter(!is.na(Gene) && !is.na(gene_length))
noncore.genes <- core.gene.assignments %>% filter(panortholog == FALSE) %>%
    inner_join(REL606.genes) %>% ## make sure these genes passed the filters
    filter(!is.na(Gene) && !is.na(gene_length))

core.gene.result <- Run.All.LTEE.Poisson.module.analysis(
    core.genes$Gene)

noncore.gene.result <- Run.All.LTEE.Poisson.module.analysis(
    noncore.genes$Gene)

adhesin.genes <- REL606.genes %>%
    filter(str_detect(product, "adhesin"))

adhesin.gene.result <- Run.All.LTEE.Poisson.module.analysis(
    adhesin.genes$Gene)

## get proteome sector assignments from Hui et al. 2015 Supplementary Table 2.
## I saved a reduced version of the data.
proteome.assignments <- read.csv(
    "../data/Hui-2015-proteome-sector-assignments.csv",
    as.is=TRUE) %>%
    inner_join(REL606.genes) %>% ## make sure these genes passed the filters
    filter(!is.na(Gene) && !is.na(gene_length))

proteome.sector.poisson.analysis.helper <- function(split.df) {
    Run.All.LTEE.Poisson.module.analysis(split.df$Gene) %>%
        ## add back the sector annotation.
        mutate(Sector.assigned = unique(split.df$Sector.assigned)) %>%
        relocate(Sector.assigned)
}

proteome.sector.results <- proteome.assignments %>%
    split(list(.$Sector.assigned)) %>%
    map_dfr(.f = proteome.sector.poisson.analysis.helper) %>%
    arrange(min.tail.fdr.corrected.pois.prob) %>%
    select(-background.mut.density, -pois.lambda)

## Get eigengene sector assignments from 
## Supplementary File 1 of Wytock and Motter (2018).
## I saved a reduced version of the data.

eigengenes <- read.csv('../data/Wytock2018-eigengenes.csv',as.is=TRUE) %>%
    inner_join(REL606.genes)

eigengene.poisson.analysis.helper <- function(split.df) {
    Run.All.LTEE.Poisson.module.analysis(split.df$Gene) %>%
        ## add back the eigengene annotation.
        mutate(Eigengene = unique(split.df$Eigengene)) %>%
        relocate(Eigengene)
}

eigengene.results <- eigengenes %>%
    split(list(.$Eigengene)) %>%
    map_dfr(.f = eigengene.poisson.analysis.helper) %>%
    arrange(min.tail.fdr.corrected.pois.prob) %>%
    select(-background.mut.density, -pois.lambda)

## Now look at genes that are regulated within Imodulons.
## I expect relaxed or purifying selection overall.

## this file was generated by reformat-I-modulons.py
genes.to.Imodulons <- read.csv("../results/gene-modules/genes-to-I-modulons.csv")
Imodulon.regulated <- genes.to.Imodulons %>% filter(!(is.na(Gene))) %>%
    ## remove any I-modulon regulators from the I-modulon genes to
    ## make an orthogonal comparison.
    filter(!(Gene %in% Imodulon.regulators$regulator))

Imodulon.regulated.poisson.analysis.helper <- function(split.df) {
    I.modulon.vec <- split.df$Gene
    Run.All.LTEE.Poisson.module.analysis(I.modulon.vec) %>%
        ## add back the I-modulon annotation.
        mutate(I.modulon = unique(split.df$I.modulon)) %>%
        relocate(I.modulon)
}

Imodulon.regulated.result <- Imodulon.regulated %>%
    split(list(.$I.modulon)) %>%
    map_dfr(.f = Imodulon.regulated.poisson.analysis.helper) %>%
    arrange(min.tail.fdr.corrected.pois.prob)

significant.Imodulon.regulated <- Imodulon.regulated.result %>%
    filter(min.tail.fdr.corrected.pois.prob < 0.01) %>%
    select(-background.mut.density, -pois.lambda)

## Note that rbs operon being under "purifying selection" is an artifact
## of parallel deletions of this region. This is excluded from the final table.
write.csv("../gene-modules/S4Table.csv", significant.Imodulon.regulated)
