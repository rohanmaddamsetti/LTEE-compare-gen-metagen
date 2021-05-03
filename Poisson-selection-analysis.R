## Poisson-selection-analysis.R by 
## Nkrumah Grant and Rohan Maddamsetti.

## TODO: make Q-Q plots to compared observed and expected distributions
## under the Poisson model. This will tell us whether or not the Poisson
## background model is a good approximation for each population.

source("metagenomics-library.R")

## get the lengths of all protein-coding genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
##Do by running:
##python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
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
    ## important: inner_join only includes genes that pass REL606 annotation filters.
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

## CRITICAL BUG CHECK TODO: make sure results don't depend on which genes are considered.
all.gene.mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    ## the full join here includes everything, including mutations in the metagenomics
    ## that don't map to the loci in REL606.genes.
    full_join(REL606.genes) %>%
    filter(Gene!='intergenic')

## CRITICAL BUG CHECK TODO: go through this list and make sure the filtered genes
## are OK (not filtered by accident)
## and make sure that these don't show up as purifying selection or anything weird.
NA.gene.mutation.data <- all.gene.mutation.data %>% filter(!complete.cases(.))

raw.metagenomic.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE)

## CRITICAL TODO: see which of these are due to errors in masking--
## see the ECB__814 prophage in particular.
missing.or.purifying <- REL606.genes %>% filter(!(Gene %in% raw.metagenomic.data$Gene))

## here is the pattern to use to check these potential bugs.
raw.metagenomic.data %>% filter(str_detect(Gene, "^ECB_008"))
raw.metagenomic.data %>% filter(str_detect(Gene, "^ECB_015"))

## CRITICAL TODO: INCLUDE THESE GENES in REL606.IDs
in.metagenomics.but.not.REL606.IDs <- raw.metagenomic.data %>%
    select(Gene) %>% distinct() %>% filter(!(Gene %in% REL606.genes$Gene))





## count number of mutated genes in each of the populations
## and then calculate the overall density of mutations
## since some regions are masked in the metagenomic data, use the total length of genes
## in REL606.genes as the normalizing factor.
TOTAL_GENE_BP <- sum(REL606.genes$gene_length)

pop.sum.mutations <- gene.mutation.data %>%
  group_by(Population) %>%
  summarise(pop.mutation.total = n(), background.mut.density = pop.mutation.total/TOTAL_GENE_BP)             

pop.poisson.data.analysis <- gene.mutation.data %>%
    group_by(Population,Gene,locus_tag,gene_length) %>%
    ## filter out synonymous mutations.
    filter(Annotation != "synonymous") %>%
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
    ## Here I hard-coded the poisson distribution as presented in the Kinnersley paper.
    ## One can also use the dpois function in base R. The results are the same.
    mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
    mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                   lower.tail=FALSE)) %>%
    ## control for the number of statistical comparisons being made:
    mutate(fdr.corrected.lower.pois.prob = p.adjust(lower.pois.prob, method = "fdr")) %>%
    mutate(fdr.corrected.upper.pois.prob = p.adjust(upper.pois.prob, method = "fdr")) %>%
    mutate(min.tail.fdr.corrected.pois.prob = pmin(fdr.corrected.lower.pois.prob,
                                                   fdr.corrected.upper.pois.prob)) %>%
    arrange(min.tail.fdr.corrected.pois.prob) %>%
    filter(min.tail.fdr.corrected.pois.prob < 0.01) %>%
    select(Population, Gene, locus_tag, gene_length, mut.count, product,
           background.mut.density,pois.lambda,fdr.corrected.lower.pois.prob,
           fdr.corrected.upper.pois.prob,min.tail.fdr.corrected.pois.prob)


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
    left_join(gene.mutation.data.total, by = character()) %>% ## add in the background mut density
    mutate(pois.lambda = background.mut.density * gene_length) %>%
    ## calculate poisson p-value:
    ## Here I hard-coded the poisson distribution as presented in the Kinnersley paper.
    ## One can also use the dpois function in base R. The results are the same.
    mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
    mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                   lower.tail=FALSE)) %>%
    ## control for the number of statistical comparisons being made:
    mutate(fdr.corrected.lower.pois.prob = p.adjust(lower.pois.prob, method = "fdr")) %>%
    mutate(fdr.corrected.upper.pois.prob = p.adjust(upper.pois.prob, method = "fdr")) %>%
    mutate(min.tail.fdr.corrected.pois.prob = pmin(fdr.corrected.lower.pois.prob,
                                                   fdr.corrected.upper.pois.prob)) %>%
    arrange(min.tail.fdr.corrected.pois.prob) %>%
    filter(min.tail.fdr.corrected.pois.prob < 0.01) %>%
    select(Gene, locus_tag, gene_length, mut.count, product,
           background.mut.density,pois.lambda,fdr.corrected.lower.pois.prob,
           fdr.corrected.upper.pois.prob,min.tail.fdr.corrected.pois.prob)

### These results look reasonable to me -- HOWEVER -- it looks like some of the purifying
### selection results (like on prophage around ECB_00845) may be artifacts of regions
### which were masked in the metagenomics data analysis, but not properly documented.
### so CAREFULLY go through these regions for debugging, using the code above!!!!!!!!



#####################

## do the Poisson analysis of modules in STIMS here.

## TODO: should I filter out synonymous mutations or not? BE MINDFUL OF HOW THIS
## AFFECTS THE RESULTS!!!
All.LTEE.Poisson.module.analysis <- function(gene.vec, gene.mutation.data, REL606.genes) {

    TOTAL_GENE_BP <- sum(REL606.genes$gene_length)
    
    ## summary over all LTEE.
    gene.mutation.data.total <- gene.mutation.data %>%
        mutate(mutation.total = nrow(.)) %>%
        mutate(background.mut.density = mutation.total/TOTAL_GENE_BP) %>%
        select(mutation.total,background.mut.density) %>% distinct()

    gene_vec_length <- sum(filter(REL606.genes,Gene %in% gene.vec)$gene_length)
    
    gene.mutation.data %>%
        ## filter for genes in gene.vec.
        filter(Gene %in% gene.vec) %>%
        ## filter out synonymous mutations. ??
##        filter(Annotation != "synonymous") %>%
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
        ## Here I hard-coded the poisson distribution as presented in the Kinnersley paper.
        ## One can also use the dpois function in base R. The results are the same.
        mutate(lower.pois.prob = ppois(q = mut.count,lambda = pois.lambda)) %>%
        mutate(upper.pois.prob = ppois(q = mut.count,lambda = pois.lambda,
                                       lower.tail=FALSE)) %>%
        ## control for the number of statistical comparisons being made:
        mutate(fdr.corrected.lower.pois.prob = p.adjust(lower.pois.prob, method = "fdr")) %>%
        mutate(fdr.corrected.upper.pois.prob = p.adjust(upper.pois.prob, method = "fdr")) %>%
        mutate(min.tail.fdr.corrected.pois.prob = pmin(fdr.corrected.lower.pois.prob,
                                                       fdr.corrected.upper.pois.prob)) %>%
        arrange(min.tail.fdr.corrected.pois.prob) %>%
        select(background.mut.density,pois.lambda,fdr.corrected.lower.pois.prob,
               fdr.corrected.upper.pois.prob,min.tail.fdr.corrected.pois.prob)   
}

## This is the function that I will use throughout, since the last two arguments
## are always the same.
Run.All.LTEE.Poisson.module.analysis <- partial(All.LTEE.Poisson.module.analysis,
                                                gene.mutation.data=gene.mutation.data,
                                                REL606.genes=REL606.genes,
                                                ... = )

neutral.genes <- read.csv("../data/neutral_compilation.csv", header=TRUE,as.is=TRUE) %>%
    ## make sure that only loci that pass filters are included in the analysis.
    filter(Gene %in% REL606.genes$Gene)

neutral.gene.result <- Run.All.LTEE.Poisson.module.analysis(neutral.genes$Gene)


## Note: this is a dataframe with one column, called regulator.
Imodulon.regulators <- read.csv("../data/rohans-I-modulons-to-regulators.csv", as.is =T) %>%
    select(regulator) %>% drop_na() %>% distinct()
    

Imodulon.regulator.result <- Run.All.LTEE.Poisson.module.analysis(Imodulon.regulators$regulator)

couce.essential.genes <- read.csv("../data/Couce2017-LTEE-essential.csv")


couce.essential.result <- Run.All.LTEE.Poisson.module.analysis(couce.essential.genes$Gene)


## get mutation parallelism in the LTEE genomes published in Tenaillon et al. (2016).
## This is used for filtering essential genes
## in the purifying selection control analysis,
## and is used in the positive selection control analysis.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv') %>%
    ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)    

hypermut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv') %>%
    ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)    

## top genes in non-mutators.
top.nonmut.genomics <- top_n(nonmut.genomics, 50, wt=G.score)
## top genes in hypermutators.
top.hypermut.genomics <- top_n(hypermut.genomics, 50, wt=G.score)

top.nonmut.result <- Run.All.LTEE.Poisson.module.analysis(top.nonmut.genomics$Gene.name)

top.hypermut.result <- Run.All.LTEE.Poisson.module.analysis(top.hypermut.genomics$Gene.name)

## core genes
core.gene.assignments <- read.csv("../data/Maddamsetti2017-core-summary.csv")
core.genes <- core.gene.assignments %>% filter(panortholog == TRUE) %>%
    left_join(REL606.genes)
noncore.genes <- core.gene.assignments %>% filter(panortholog == FALSE) %>%
    left_join(REL606.genes)

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
proteome.assignments <- read.csv('../data/Hui-2015-proteome-sector-assignments.csv',
                                 as.is=TRUE) %>%
    inner_join(REL606.genes)
## add proteome assignment to gene mutation.data.
sector.mut.data <- inner_join(gene.mutation.data, proteome.assignments)
##six sectors:  "A" "S" "O" "U" "R" "C"
A.sector <- filter(proteome.assignments,Sector.assigned=='A')
S.sector <- filter(proteome.assignments,Sector.assigned=='S')
O.sector <- filter(proteome.assignments,Sector.assigned=='O')
U.sector <- filter(proteome.assignments,Sector.assigned=='U')
R.sector <- filter(proteome.assignments,Sector.assigned=='R')
C.sector <- filter(proteome.assignments,Sector.assigned=='C')

A.sector.result <- All.LTEE.Poisson.module.analysis(
    A.sector$Gene,
    gene.mutation.data,
    REL606.genes)

S.sector.result <- All.LTEE.Poisson.module.analysis(
    S.sector$Gene,
    gene.mutation.data,
    REL606.genes)

O.sector.result <- All.LTEE.Poisson.module.analysis(
    O.sector$Gene,
    gene.mutation.data,
    REL606.genes)

U.sector.result <- All.LTEE.Poisson.module.analysis(
    U.sector$Gene,
    gene.mutation.data,
    REL606.genes)

R.sector.result <- All.LTEE.Poisson.module.analysis(
    R.sector$Gene,
    gene.mutation.data,
    REL606.genes)

C.sector.result <- All.LTEE.Poisson.module.analysis(
    C.sector$Gene,
    gene.mutation.data,
    REL606.genes)

## get eigengene sector assignments from Wytock and Motter (2018) Supplementary File 1.
## I saved a reduced version of the data.

eigengenes <- read.csv('../data/Wytock2018-eigengenes.csv',as.is=TRUE) %>%
    inner_join(REL606.genes)

eigengene1 <- filter(eigengenes, Eigengene==1)
eigengene2 <- filter(eigengenes, Eigengene==2)
eigengene3 <- filter(eigengenes, Eigengene==3)
eigengene4 <- filter(eigengenes, Eigengene==4)
eigengene5 <- filter(eigengenes, Eigengene==5)
eigengene6 <- filter(eigengenes, Eigengene==6)
eigengene7 <- filter(eigengenes, Eigengene==7)
eigengene8 <- filter(eigengenes, Eigengene==8)
eigengene9 <- filter(eigengenes, Eigengene==9)

eigengene1.result <- All.LTEE.Poisson.module.analysis(
    eigengene1$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene2.result <- All.LTEE.Poisson.module.analysis(
    eigengene2$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene3.result <- All.LTEE.Poisson.module.analysis(
    eigengene3$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene4.result <- All.LTEE.Poisson.module.analysis(
    eigengene4$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene5.result <- All.LTEE.Poisson.module.analysis(
    eigengene5$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene6.result <- All.LTEE.Poisson.module.analysis(
    eigengene6$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene7.result <- All.LTEE.Poisson.module.analysis(
    eigengene7$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene8.result <- All.LTEE.Poisson.module.analysis(
    eigengene8$Gene,
    gene.mutation.data,
    REL606.genes)

eigengene9.result <- All.LTEE.Poisson.module.analysis(
    eigengene9$Gene,
    gene.mutation.data,
    REL606.genes)

## Now look at genes that are regulated within Imodulons.
## I expect relaxed or purifying selection overall.

## this file was generated by reformat-I-modulons.py
genes.to.Imodulons <- read.csv("../results/gene-modules/genes-to-I-modulons.csv")
Imodulon.regulated <- genes.to.Imodulons %>% filter(!(is.na(Gene)))

Run.All.LTEE.Poisson.analysis.on.Imodulons <- function(Imodulon.regulated) {

    helper.func <- function(df) {
        I.modulon.vec <- df$Gene
        return(Run.All.LTEE.Poisson.module.analysis(I.modulon.vec))
    }
    
    Imodulon.regulated %>%
        split(list(.$I.modulon)) %>%
        map_dfr(.f=helper.func)
}

Imodulon.regulated.result <- Run.All.LTEE.Poisson.analysis.on.Imodulons(Imodulon.regulated)

#################################

## Rohan messing around with results here.


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

## example: rhsB, rhsD, rhsE regions are hypermutable due to gene conversion I believe.

## POTENTIAL TODO: more sophisticated/complicated background mutation density for
## poisson model. Could fit a local mutation density by binning the genome, like
## in divergent mutation bias paper, and some of the STIMS code allowing for
## regional mutation variation, etc. May not be worth the extra complexity though.
