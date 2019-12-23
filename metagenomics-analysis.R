## metagenomics-analysis.R by Rohan Maddamsetti.

## Basic premise.
## count the cumulative number of stars over time, and plot.
## examine different kinds of mutations and genes.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## TODO: infer cohorts using Haixu Tang's new algorithm.
## Then ask whether cohorts show functional enrichment using
## STRING annotation.

## TODO: It could be quite interesting to examine genes in prophage that
## Jeff Barrick has annotated as misc_regions.

## TODO: measure the concentration/density of mutation as an empirical CDF.

library(tidyverse)
library(DescTools) ## for G-test implementation in the multinomial selection test.

## get the lengths of all genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
##Do by running:
##python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
mutate(gene_length=strtoi(gene_length))

## IMPORTANT BUG: SEEMS LIKE BEN MASKED SOME REGIONS--
## including prophage starting at ECB_00814 --
## THAT I MISS!!!
## I can double-check my purifying selection results against the
## LTEE genomics data at barricklab.org/shiny/LTEE-Ecoli
## to verify that those genes indeed don't have any mutations.


## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000)

gene.mutation.data <- inner_join(mutation.data,REL606.genes)

## It turns out that some gene names map to multiple genes!!!
duplicate.genes <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    group_by(Gene) %>%
    summarize(checkme=length(unique(gene_length))) %>%
    filter(checkme>1)

## filter those duplicates.
gene.mutation.data <- filter(gene.mutation.data,
                             !(Gene %in% duplicate.genes$Gene))

## write out gene sets to examine using the STRING database.
write.csv(gene.mutation.data, "../results/gene-LTEE-metagenomic-mutations.csv")


########################################################################
## investigate dS across the genome in the metagenomics data.
## revamp code from my 2015 Mol. Biol. Evol. paper.
## in short, cannot reject null model that dS is uniform over the genome.
## and dN fits the null extremely well-- even better than dS!
## Probably because there are 3 times as many dN as dS throughout the
## experiment.
ks.analysis <- function (the.data,REL606.genes) {
  ## For each set of data (all data, non-mutators, MMR mutators, mutT mutators)
  ## do the following: 1) make a uniform cdf on mutation rate per base.
  ## 2) make an empirical cdf of mutations per gene.
  ## do K-S tests for goodness of fit of the empirical cdf with the cdfs for
  ## the uniform cdf and thetaS cdf hypotheses.

  hit.genes.df <- the.data %>%
  group_by(locus_tag,Gene,gene_length) %>%
  summarize(hits=n()) %>%
  ungroup()

  ## have to do it this way, so that zeros are included.
  hit.genes.df <- full_join(REL606.genes,hit.genes.df) %>% replace_na(list(hits=0)) %>%
  arrange(desc(gene_length))
    
  ## Calculate the empirical distribution of synonymous substitutions per gene.
  hit.genes.length <- sum(hit.genes.df$gene_length)
  mutation.total <- sum(hit.genes.df$hits)
  empirical.cdf <- cumsum(hit.genes.df$hits)/mutation.total
  ## Null hypothesis: probability of a mutation per base is uniform.
  null.cdf <- cumsum(hit.genes.df$gene_length)/hit.genes.length

  ## Do Kolmogorov-Smirnov tests for goodness of fit.
  print(ks.test(empirical.cdf, null.cdf, simulate.p.value=TRUE))

    results.to.plot <- data.frame(locus_tag=hit.genes.df$locus_tag,
                                  Gene=hit.genes.df$Gene,
                                  gene_length=hit.genes.df$gene_length,
                                  empirical=empirical.cdf,
                                  null=null.cdf)
  return(results.to.plot)
}

make.KS.Figure <- function(the.results.to.plot) { 
  ## for plotting convenience, add an index to the data frame.
    the.results.to.plot$index <- 1:length(the.results.to.plot$locus_tag)
  
    plot <- ggplot(the.results.to.plot, aes(x=index)) +
        geom_line(aes(y=empirical), colour="red") + 
        geom_line(aes(y=null), linetype=2) + 
        scale_x_continuous('Genes ranked by length',limits=c(0,4400)) +
        scale_y_continuous('Cumulative proportion of mutations',limits=c(0,1)) +
        theme_classic() +
        theme(axis.title=element_text(size=18),axis.text=element_text(size=12))
    
    return(plot)
}

## examine all mutations over the genome.
cumsum.all.over.metagenome <- ks.analysis(gene.mutation.data,REL606.genes)
make.KS.Figure(cumsum.all.over.metagenome)

## examine dS over the genome.
gene.dS.mutation.data <- gene.mutation.data %>%
filter(Annotation=='synonymous')

cumsum.dS.over.metagenome <- ks.analysis(gene.dS.mutation.data,REL606.genes)
make.KS.Figure(cumsum.dS.over.metagenome)

## examing dN over the genome.
gene.dN.mutation.data <- gene.mutation.data %>%
filter(Annotation=='missense')

## dN fits the null even better than dS! But note that
## there are 18493 dN in the data, and 6792 dS in the data.
## so the better fit is probably best explained by the larger sample size.
cumsum.dN.over.metagenome <- ks.analysis(gene.dN.mutation.data,REL606.genes)
make.KS.Figure(cumsum.dN.over.metagenome)

## let's look at nonsense mutations.
gene.nonsense.mutation.data <- gene.mutation.data %>%
filter(Annotation=='nonsense')

## nonsense mutations don't fit the null expectation.
## opposite trend of indels or IS elements, though!
## not sure why.
cumsum.nonsense.over.metagenome <- ks.analysis(gene.nonsense.mutation.data,REL606.genes)
make.KS.Figure(cumsum.nonsense.over.metagenome)

## now let's look at all mutations except for dS.
gene.except.dS.mutation.data <- gene.mutation.data %>%
filter(Gene!='intergenic') %>%
filter(Annotation!='synonymous')

cumsum.no.dS.over.metagenome <- ks.analysis(gene.except.dS.mutation.data,REL606.genes)
make.KS.Figure(cumsum.no.dS.over.metagenome)

## There's a significant difference between dS and non-dS
## distribution over the genome:
## p = 0.001. I feel I found this by fishing... but still
## significant by Bonferroni correcting by the different graphs I made
## in this section (9 tests)
ks.test(cumsum.dS.over.metagenome$empirical,
        cumsum.no.dS.over.metagenome$empirical,
        simulate.p.value=TRUE)
## Two hypotheses for these results:
## H1: mutation hotspots for indels and IS elements.
## H2: purifying selection against indels and IS elements.

## but this test is not significant in comparing dS to dN.
## so driven by distribution of non-point mutations over the genome,
## since I excluded intergenic mutations.
ks.test(cumsum.dS.over.metagenome$empirical,
        cumsum.dN.over.metagenome$empirical,
        simulate.p.value=TRUE)

## examine structural mutations (IS elements) affecting genes.
gene.sv.mutation.data <- gene.mutation.data %>%
filter(Gene!='intergenic') %>%
filter(Annotation=='sv')

## RESULT: longer genes are depleted in IS insertions!!
cumsum.sv.over.metagenome <- ks.analysis(gene.sv.mutation.data,REL606.genes)
make.KS.Figure(cumsum.sv.over.metagenome)

## examine indels.
gene.indel.mutation.data <- gene.mutation.data %>%
filter(Gene!='intergenic') %>%
filter(Annotation=='indel')

## RESULT: longer genes are depleted in indels!
cumsum.indel.over.metagenome <- ks.analysis(gene.indel.mutation.data,REL606.genes)
make.KS.Figure(cumsum.indel.over.metagenome)

## Therefore, structural mutations and indels are probably what are driving the
## differences that I saw in aerobic and anerobic mutations before.

## let's combine IS (structural mutations), indels, and nonsense mutations,
## as a proxy for purifying selection.
gene.nonsense.sv.indels.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation %in% c('indel', 'sv', 'nonsense'))

cumsum.nonsense.sv.indels.over.metagenome <- ks.analysis(gene.nonsense.sv.indels.mutation.data, REL606.genes)
make.KS.Figure(cumsum.nonsense.sv.indels.over.metagenome)

## TODO: I WILL HAVE TO DO MORE CAREFUL WORK TO EXAMINE THE EFFECTS OF MUTATION BIAS.
## look at structural variation. can I distinguish between
## mutation hotspots vs. selection?
## Right now, I don't think I can.

## TODO: look at mutation density over the chromosome. Any wave patterns,
## as in recent work from Pat Foster's group? Per my reading of the
## abstract of their most recent paper, perhaps the uniformity
## of dS over LTEE genomes reflect unwinding/loosening of chromosomal
## proteins packing up the DNA. Just speculation.

## TODO: I need to examine mutation biases in further detail, in order to
## get past peer review.

## cross-check these results with 50K genomes.
## These were downloaded from Jeff Barrick's shiny app:
## https://barricklab.org/shiny/LTEE-Ecoli.

gene.50K.mutations <- read.csv("../data/Gen50000_allMutations.csv",
                              header=TRUE,
                              as.is=TRUE) %>%
    filter(clone=='A') %>%
    rename(Gene=gene_name) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

## cross-check with dS in 50K genomes.
dS.50K <- filter(gene.50K.mutations, snp_type=='synonymous')
cumsum.dS.genome.50K <- ks.analysis(dS.50K, REL606.genes)
make.KS.Figure(cumsum.dS.genome.50K)

## cross-check with everything except dS in 50K genomes.
no.dS.50K <- filter(gene.50K.mutations, snp_type!='synonymous')
cumsum.no.dS.genome.50K <- ks.analysis(no.dS.50K, REL606.genes)
make.KS.Figure(cumsum.no.dS.genome.50K)

## cross-check with dN in 50K genomes.
dN.50K <- filter(gene.50K.mutations, snp_type!='nonsynonymous')
cumsum.dN.genome.50K <- ks.analysis(dN.50K, REL606.genes)
make.KS.Figure(cumsum.dN.genome.50K)

## cross-check IS elements in the 50K genomes.
IS.50K <- filter(gene.50K.mutations, mutation_category=="mobile_element_insertion")
cumsum.IS.genome.50K <- ks.analysis(IS.50K,REL606.genes)
make.KS.Figure(cumsum.IS.genome.50K)

## cross-check indels in the 50K genomes.
indels.50K <- filter(gene.50K.mutations, mutation_category=="small_indel")
cumsum.indels.genome.50K <- ks.analysis(indels.50K,REL606.genes)
make.KS.Figure(cumsum.indels.genome.50K)

## cross-check nonsense SNPs in the 50K genomes.
## OPPOSITE trend!!!
nonsense.50K <- filter(gene.50K.mutations, mutation_category=="snp_nonsense")
cumsum.nonsense.genome.50K <- ks.analysis(nonsense.50K,REL606.genes)
make.KS.Figure(cumsum.nonsense.genome.50K)

## large deletions in 50K genomes.
largedeletions.50K <- filter(gene.50K.mutations, mutation_category=="large_deletion")
cumsum.largedeletions.genome.50K <- ks.analysis(largedeletions.50K,REL606.genes)
make.KS.Figure(cumsum.largedeletions.genome.50K)

## The dynamics of mutations under relaxed selection should be
## dragged by mutations under positive selection.
## Also see Ville Mustonen's paper on 'emergent neutrality'.
## so a better null hypothesis: dynamics of ALL mutations are driven by
## positive selection, either as hitchhikers or as drivers.

## Does the Ornstein-Uhlenbeck process require that effects of each
## mutation is uncorrelated? If mutations affecting anaerobic fitness
## are uncorrelated with mutations affecting aerobic fitness,
## then would this Genotype-Phenotype Map be inconsistent
## with Nkrumah's results?

#################################################################################

## For every gene, calculate my multinomial test for fit to a "neutral" model of
## hits per gene, based on the number of mutations in each population per mutation class.
## combine this in a dataframe.

## Then,
## Calculate the total density of mutations for each gene
## (for different classes of mutations),
## and combine in a dataframe.

## I want to use both sets of these statistics, on a per-gene basis, in order to
## see if any of these features relate to the rates at which various classes of
## mutations occur in different sets of genes.


################
## calculate the probability of a given configuration of parallel mutations
## in the metagenomes across populations.
## This seems to be a neat test for positive selection!
## Seems to give the same answer as results in Tenaillon Nature paper.
## TODO: compare to Supplementary Table S3 in Good et al. Nature paper.
gene.multinom.probability <- function (mutation.data,gene,mut.class.vec) {
  population.probs <- mutation.data %>% group_by(Population) %>%
  summarize(total.muts=n()) %>% mutate(prob=total.muts/sum(total.muts))

  gene.data <- filter(mutation.data,Gene==gene,Annotation %in% mut.class.vec)
    
  gene.data.counts <- gene.data %>% group_by(Population) %>%
  summarize(muts=n())

  if(nrow(gene.data.counts) == 0) return(NA)
  
  ## hack to get lengths of vectors matching for the GTest.
  full.pop.column <- data.frame(Population=population.probs$Population,
                                stringsAsFactors=FALSE)
  stat.df <- full_join(full.pop.column,gene.data.counts,by='Population')
  stat.df[is.na(stat.df)] <- 0

  ## Use a G-test from the DescTools package because
  ## an exact multinomial test is too slow.
  GTest(x=stat.df$muts, p=population.probs$prob)$p.value
}

## draw every gene in the genome. calculate log probability and rank.
genes.vec <- unique(REL606.genes$Gene)
all.muts.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("missense", "synonymous", "sv", "indel", "nonsense"))})
dN.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("missense"))})
dS.pval.vec <- sapply(genes.vec,function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("synonymous"))})
nonsense.sv.indel.pval.vec <- sapply(genes.vec, function(gene) {gene.multinom.probability(gene.mutation.data,gene,c("sv", "indel", "nonsense"))})

## IMPORTANT: KEEP NA VALUES.
multinom.result <- data.frame(Gene=genes.vec,
                              dN.pvalue=dN.pval.vec,
                              dS.pvalue=dS.pval.vec,
                              all.mut.pvalue=all.muts.pval.vec,
                              nonsense.sv.indel.pvalue=nonsense.sv.indel.pval.vec) %>%
    arrange(dN.pvalue) %>%
    mutate(index=1:n()) %>%
    mutate(dS.qvalue=p.adjust(dS.pvalue,'fdr')) %>%
    mutate(all.mut.qvalue=p.adjust(all.mut.pvalue,'fdr')) %>%
    mutate(nonsense.sv.indel.qvalue=p.adjust(nonsense.sv.indel.pvalue,'fdr')) %>%
    mutate(dN.qvalue=p.adjust(dN.pvalue,'fdr')) %>%
    ## keep the annotation in REL606.genes.
    full_join(REL606.genes)
### IMPORTANT BUG! GENES THAT SHOULD HAVE BEEN MASKED, WITH NO
### MUTATION CALLS IN THE GOOD ET AL. DATA, REAPPEAR IN REL606 GENES!!!
## PROPHAGE STARTING AT ECB_00814 is a good place to start debugging.


## write multinomial selection test results to file.
write.csv(multinom.result,file='../results/multinomial_test_for_selection.csv')

## take a look at these distributions. The top genes are all under strong
## positive selection in the LTEE, as reported by Tenaillon et al.
multinom.plot <- ggplot(multinom.result,aes(x=index,y=-log(dN.pvalue),label=Gene,color)) + geom_point() + theme_classic()

multinom.plot2 <- ggplot(arrange(multinom.result,all.mut.pvalue),aes(x=index,y=-log(all.mut.pvalue),label=Gene,color)) + geom_point() + theme_classic()

## Now, let's examine the results of this test.

## interesting! cpsG is significant in terms of dS!
sig.multinom.dS <- filter(multinom.result,dS.pvalue<0.05)

## most of these genes are already known. Some aren't.
sig.multinom.dN <- filter(multinom.result,dN.pvalue<0.05)

sig.multinom.all.muts <- filter(multinom.result,all.mut.pvalue<0.05)
sig.multinom..nonsense.sv.indels <- filter(multinom.result,nonsense.sv.indel.pvalue<0.05)

## look at genes which are significant based on this test, but
## not based on parallelism alone, from Ben Good Table S3.
## these may be candidates for contingency.
Good.significant.genes <- read.csv('../ltee-metagenomics-paper/nature24287-s5.csv.txt')
contingency.candidates <- filter(sig.multinom.dN,!(Gene %in% Good.significant.genes$Gene))
## of these genes, hflB, topA are in the Tenaillon top G-scoring
## genes. The rest are marginally significant.
## ECB_03460 is interesting though-- this is an adhesin! Why is this here?

## worth querying these genes (and the significant dS gene!)
## in the Good metagenomic data to see if there's anything
## interesting.

################
## Examine the distribution of various classes of mutations across genes in the
## genomics or metagenomics data. Which genes are enriched? Which genes are depleted?
## Then, can look at the annotation of these genes in STRING.
calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    density.df <- gene.mutation.data %>%
        filter(Annotation %in% mut_type_vec) %>%
        filter(Gene!= "intergenic") %>%
        mutate(Gene=as.factor(Gene)) %>%
        group_by(Gene,gene_length) %>%
        summarize(mut.count=n()) %>%
        ungroup() %>%
        mutate(density=mut.count/gene_length) %>%
        arrange(desc(density))
    return(density.df)
}

## just rank gene.mutation.density to see what it looks like!
## this is already a good indication of selection,
## and this is pretty much what Ben Good does in his Table S3.

all.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data,
    c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense")) %>%
    rename(all.mut.count = mut.count) %>%
    rename(all.mut.density = density)

dN.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data, c("missense")) %>%
    rename(dN.mut.count = mut.count) %>%
    rename(dN.mut.density = density)

dS.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data, c("synonymous")) %>%
    rename(dS.mut.count = mut.count) %>%
    rename(dS.mut.density = density)

sv.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data,c("sv")) %>%
    rename(sv.mut.count = mut.count) %>%
    rename(sv.mut.density = density)

indel.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data,c("indel")) %>%
    rename(indel.mut.count = mut.count) %>%
    rename(indel.mut.density = density)

nonsense.indel.sv.density <- calc.gene.mutation.density(
    gene.mutation.data,c("sv", "indel", "nonsense")) %>%
    rename(nonsense.indel.sv.mut.count = mut.count) %>%
    rename(nonsense.indel.sv.mut.density = density)

all.except.dS.density <- calc.gene.mutation.density(
    gene.mutation.data,c("sv", "indel", "nonsense", "missense")) %>%
    rename(all.except.dS.mut.count = mut.count) %>%
    rename(all.except.dS.mut.density = density)

## combine these into one dataframe.
gene.mutation.densities <- REL606.genes %>%
    full_join(all.mutation.density) %>%
    full_join(dN.mutation.density) %>%
    full_join(dS.mutation.density) %>%
    full_join(sv.mutation.density) %>%
    full_join(indel.mutation.density) %>%
    full_join(nonsense.indel.sv.density) %>%
    full_join(all.except.dS.density)

## replace NAs with zeros.
gene.mutation.densities[is.na(gene.mutation.densities)] <- 0
gene.mutation.densities <- tbl_df(gene.mutation.densities)

## how does the multinomial rest compare to the number of mutations in each gene?
## Control for gene length, and see if there's a correlation between log(p)
## and density of mutations in each gene.
## there's a super strong relationship with dN p-value: p < 10^-15.
multinom.result.and.density <- inner_join(gene.mutation.densities,multinom.result)
cor.test(multinom.result.and.density$dN.mut.density,multinom.result.and.density$dN.pvalue)
## Also strong correlation with dS p-value: p < 10^-14.
cor.test(multinom.result.and.density$dS.mut.density,multinom.result.and.density$dS.pvalue)

#################################################
## PURIFYING SELECTION RESULTS.
## VERY PROMISING BUT WILL NEED TO RIGOROUSLY CHECK.

## Let's look at distribution of SV and indels across genes in the
## metagenomics data. Which genes are enriched? Which genes are depleted?
## Then, look at the annotation of these genes in STRING.

## get genes with no sv.
no.sv.genes <- filter(gene.mutation.densities, sv.mut.count == 0) %>%
    arrange(desc(gene_length))
## get genes with no indels.
no.indel.genes <- filter(gene.mutation.densities, indel.mut.count == 0) %>%
    arrange(desc(gene_length))

## IMPORTANT: EXAMINE GENES WITH ZERO INDEL, SV, NONSENSE MUTATIONS.
## THESE ARE THE MOST DEPLETED.
## THIS LOOKS REAL!!!
no.nonsense.indel.sv.genes <- filter(gene.mutation.densities,
                                     nonsense.indel.sv.mut.count == 0) %>%
    arrange(desc(gene_length))
write.csv(no.nonsense.indel.sv.genes,file="../results/no-nonsense-indel-sv-genes.csv")

## Try again, now disallowing missense mutations.
only.dS.allowed.genes <- filter(gene.mutation.densities,
                                all.except.dS.mut.count == 0) %>%
    arrange(desc(gene_length))
write.csv(only.dS.allowed.genes,file="../results/only-dS-allowed-genes.csv")


## As an added confirmation, let's compare essentiality from KEIO collection to these
## gene sets.
KEIO.data <- read.csv("../data/KEIO_Essentiality.csv", header=TRUE,as.is=TRUE) %>%
    select(-JW_id)

KEIO.gene.mutation.densities <- left_join(gene.mutation.densities,KEIO.data)


purifying1 <- KEIO.gene.mutation.densities %>% mutate(maybe.purifying=ifelse(nonsense.indel.sv.mut.density==0,TRUE,FALSE))

purifying1.plot <- ggplot(purifying1,aes(x=maybe.purifying,y=Score)) + theme_classic() + geom_boxplot()
purifying1.plot

purifying2 <- purifying1 %>% filter(gene_length > 1000)

purifying2.plot <- ggplot(purifying2,aes(x=maybe.purifying,y=Score)) + theme_classic() + geom_boxplot()
purifying2.plot

purifying3 <- KEIO.gene.mutation.densities %>% mutate(maybe.purifying=ifelse(all.except.dS.mut.density==0,TRUE,FALSE))

purifying3.plot <- ggplot(purifying3,aes(x=maybe.purifying,y=Score)) + theme_classic() + geom_boxplot()
purifying3.plot

## potentially purifying genes have a higher KEIO essentially score, as we would hope.
wilcox.test(x=filter(purifying1,maybe.purifying==TRUE)$Score,filter(purifying1,maybe.purifying==FALSE)$Score)

wilcox.test(x=filter(purifying2,maybe.purifying==TRUE)$Score,filter(purifying2,maybe.purifying==FALSE)$Score)

wilcox.test(x=filter(purifying3,maybe.purifying==TRUE)$Score,filter(purifying3,maybe.purifying==FALSE)$Score)


pur1 <- filter(purifying1,maybe.purifying==TRUE) ## looks like a lot of false positives
pur2 <- filter(purifying2,maybe.purifying==TRUE) ## looks like a lot of false positives

## ribosomal genes... and lots of prophage? Are these addictive genes?

## See "Pervasive domestication of defective prophages by bacteria"
## by Louis-Marie Bobay for insights, or potential analyses and checks, perhaps?

## Also see: Bacterial ‘Grounded’ Prophages: Hotspots for Genetic Renovation and Innovation
## by Ramisetty and Sudhakari.

pur3 <- filter(purifying3,maybe.purifying==TRUE) %>% arrange(desc(gene_length))

## by allowing dS, I might be able to filter out false positives caused by
## a lack of reads mapping uniquely to that locus,
## such as repeats or transposases or something.
pur3.with.dS <- filter(pur3,dS.mut.count>0)

## I'll need to dig deeper into these results.


## IMPORTANT BUG: SEEMS LIKE BEN MASKED SOME REGIONS--
## including prophage starting at ECB_00814 --
## THAT I MISS!!!
## I can double-check my purifying selection results against the
## LTEE genomics data at barricklab.org/shiny/LTEE-Ecoli
## to verify that those genes indeed don't have any mutations.
##############################################
################################################################################
## Now, look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.

## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length)

cumsum.per.pop.helper.func <- function(df) {
  df %>%
  arrange(t0) %>%
  group_by(Population,Generation) %>%
  summarize(count=n()) %>%
  mutate(cs=cumsum(count))
}

calc.cumulative.muts <- function(data, normalization.constant) {
  data %>%
  split(.$Population) %>%
  map_dfr(.f=cumsum.per.pop.helper.func) %>%
  mutate(normalized.cs=cs/normalization.constant)
}

###############################################################
## Examine the rate of cumulative accumulation of mutations over time in gene sets
## selected based on the multinomial measures.

## break REL606 genes into deciles, based both on mutation density, as well as
## based on p-values from the multinomial test.

## decile = 1 means lowest density. decile = 10 means highest density.
gene.deciles.by.density <- gene.mutation.densities %>%
    mutate(all.mut.decile = ntile(all.mut.density, 10)) %>%
    mutate(dN.decile = ntile(dN.mut.density, 10)) %>%
    mutate(nonsense.indel.sv.decile = ntile(nonsense.indel.sv.mut.density, 10))

## examine density tails by dN.
median.by.dN.density.genes <- filter(gene.deciles.by.density,dN.decile %in% c(5,6))
median.by.dN.mutation.data <- filter(gene.mutation.data,Gene %in% median.by.dN.density.genes$Gene)
median.by.dN.gene.length <- sum(median.by.dN.mutation.data$gene_length)

top.by.dN.density.genes <- filter(gene.deciles.by.density,dN.decile == 10)
top.by.dN.mutation.data <- filter(gene.mutation.data,Gene %in% top.by.dN.density.genes$Gene)
top.by.dN.gene.length <- sum(top.by.dN.mutation.data$gene_length)

bottom.by.dN.density.genes <- filter(gene.deciles.by.density,dN.decile == 1)
bottom.by.dN.mutation.data <- filter(gene.mutation.data,Gene %in% bottom.by.dN.density.genes$Gene)
bottom.by.dN.gene.length <- sum(bottom.by.dN.mutation.data$gene_length)

## examine density tails by nonsense, indels, and sv.
median.by.nonindelsv.density.genes <- filter(gene.deciles.by.density,nonsense.indel.sv.decile %in% c(5,6))
median.by.nonindelsv.mutation.data <- filter(gene.mutation.data,Gene %in% median.by.nonindelsv.density.genes$Gene)
median.by.nonindelsv.gene.length <- sum(median.by.nonindelsv.mutation.data$gene_length)

top.by.nonindelsv.density.genes <- filter(gene.deciles.by.density,nonsense.indel.sv.decile == 10)
top.by.nonindelsv.mutation.data <- filter(gene.mutation.data,Gene %in% top.by.nonindelsv.density.genes$Gene)
top.by.nonindelsv.gene.length <- sum(top.by.nonindelsv.mutation.data$gene_length)

bottom.by.nonindelsv.density.genes <- filter(gene.deciles.by.density,nonsense.indel.sv.decile == 1)
bottom.by.nonindelsv.mutation.data <- filter(gene.mutation.data,Gene %in% bottom.by.nonindelsv.density.genes$Gene)
bottom.by.nonindelsv.gene.length <- sum(bottom.by.nonindelsv.mutation.data$gene_length)

## NOTE: The multinom analysis here is based ONLY on distribution of dN.

## NOTE: rate of accumulation depends on the size of the gene set,
## when looking at the left and right tails.
## (I can't remember what this comment means-- take a look at the data to figure out.)

## also compare genes in the tails of the multinomial hit distribution.
## left tail is strong selection. right tail fits the null well.

## ODD! way more mutations in top 100 or top 200 genes.
## but way fewer in the top 50! Take a look at 1-50, and 51-100 genes.

multinom.1to50.genes <- multinom.result$Gene[1:50]
multinom.50to100.genes <- multinom.result$Gene[51:100]

multinom.1to100.genes <- multinom.result$Gene[1:100]
multinom.1to200.genes <- multinom.result$Gene[1:200]


multinom.1to50.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to50.genes)

multinom.50to100.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.50to100.genes)

multinom.1to100.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to100.genes)

multinom.1to200.mutation.data <- filter(mutation.data,
                                  Gene %in% multinom.1to200.genes)

## for the bottom, get the number of genes in the multinom.result
## (for now, including NA values).

num.multinom.genes <- length(multinom.result$Gene)

multinom.bottom50to1.genes <- multinom.result$Gene[num.multinom.genes-50:num.multinom.genes]
multinom.bottom100to50.genes <- multinom.result$Gene[(num.multinom.genes-100):(num.multinom.genes-50)]
multinom.bottom100to1.genes <- multinom.result$Gene[(num.multinom.genes-100):num.multinom.genes]
multinom.bottom200to1.genes <- multinom.result$Gene[num.multinom.genes-200:num.multinom.genes]

multinom.bottom200to1.mutation.data <- filter(mutation.data,
                                              Gene %in% multinom.bottom200to1.genes)

## length of tails of multinomial distibution.
multinom.1to50.length <- sum(filter(REL606.genes, Gene %in% multinom.1to50.genes)$gene_length, na.rm=TRUE)
multinom.50to100.length <- sum(filter(REL606.genes, Gene %in% multinom.50to100.genes)$gene_length, na.rm=TRUE)
multinom.1to100.length <- sum(filter(REL606.genes, Gene %in% multinom.1to100.genes)$gene_length, na.rm=TRUE)
multinom.1to200.length <- sum(filter(REL606.genes, Gene %in% multinom.1to200.genes)$gene_length, na.rm=TRUE)

multinom.bottom50to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom50to1.genes)$gene_length, na.rm=TRUE)
multinom.bottom100to50.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom100to50.genes)$gene_length, na.rm=TRUE)
multinom.bottom100to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom100to1.genes)$gene_length, na.rm=TRUE)
multinom.bottom200to1.length <- sum(filter(REL606.genes, Gene %in% multinom.bottom200to1.genes)$gene_length, na.rm=TRUE)


#########################################################################
## Finally -- after all this prep work--
## calculate cumulative numbers of mutations in each category.
## Then, make plots to see if any interesting patterns emerge.

plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log(normalized.cs)))
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs))
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='fixed') +
        ylab('Cumulative number of mutations, normalized by target size') +
        xlab('Generations (x 10,000)')
    return(p)
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=TRUE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log(normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, size=0.2, color=my.color)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, size=0.2, color=my.color)
        }
    return(p)
}

#########################################################################

## Normalization constant calculations.
## IMPORTANT-- THIS STEP IS REALLY IMPORTANT.
## NOT CLEAR HOW TO NORMALIZE IN ORDER TO COMPARE DIFFERENT CLASSES OF MUTATIONS
## "APPLES TO APPLES".

## TODO: check difference between normalizing by gene length and normalizing
## by synonymous/nonsynonymous opportunities. The end result should be very similar.
## Then consider refactoring to just use REL606_IDs.csv to get gene lengths.

## numbers gotten by running measureTargetSize.py.
## Use these to normalize cumulative mutations over time.
## IMPORTANT TODO: THESE NUMBERS ARE ALL OFF NOW-- HAVE TO TAKE MASKING INTO ACCOUNT!!!
target.size.numbers <- read.csv('../results/target_size.csv',header=TRUE,as.is=TRUE)

total.length <- filter(target.size.numbers,set=='genome')$total_gene_length
total.synon.sites <- filter(target.size.numbers,set=='genome')$synon_sites
total.nonsynon.sites <- filter(target.size.numbers,set=='genome')$non_synon_sites

######################NOTE: DOUBLE CHECK CONSISTENT WITH measureTargetSize.py. output!!!!!
## IMPORTANT TODO: THESE NUMBERS ARE ALL OFF NOW-- HAVE TO TAKE MASKING INTO ACCOUNT!!!
####### Constants. REVAMP THIS CODE!
## from measureIntergenicTargetSize.py:
## Length of intergenic regions: 487863
intergenic.length <- 487863


########################
## let's first look at the accumulation of all sv, indels, and nonsense mutations
## in each population.

sv.indel.nonsense.mutation.data <- mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))

c.sv.indel.nonsense.mutations <- calc.cumulative.muts(sv.indel.nonsense.mutation.data,total.length)

c.sv.indel.nonsense.plot <- plot.cumulative.muts(c.sv.indel.nonsense.mutations)
c.sv.indel.nonsense.plot <- plot.cumulative.muts(c.sv.indel.nonsense.mutations,logscale=FALSE)
c.sv.indel.nonsense.plot

## Now filter for nonsense, sv, and indels in genes (of course nonsense are always in genes.)
sv.indel.nonsense.gene.mutation.data <- gene.mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))
c.sv.indel.nonsense.gene.mutations <- calc.cumulative.muts(sv.indel.nonsense.gene.mutation.data,total.length)

c.sv.indel.nonsense.gene.plot <- plot.cumulative.muts(c.sv.indel.nonsense.gene.mutations)
c.sv.indel.nonsense.gene.plot <- plot.cumulative.muts(c.sv.indel.nonsense.gene.mutations,logscale=FALSE)

c.sv.indel.nonsense.gene.plot


########################
## Examine ALL mutations, including intergenic regions.
c.mutations <- calc.cumulative.muts(mutation.data,total.length)

c.mutation.plot <- plot.cumulative.muts(c.mutations, logscale=TRUE)
c.mutation.plot2 <- plot.cumulative.muts(c.mutations, logscale=FALSE)

## Examine dN mutations.
##c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
dN.normalization.const <- total.length * total.nonsynon.sites/(total.synon.sites+total.nonsynon.sites)
c.dN.mutations <- calc.cumulative.muts(gene.dN.mutation.data,dN.normalization.const)

## Examine dS mutations.
##c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)
dS.normalization.const <- total.length * total.synon.sites/(total.synon.sites+total.nonsynon.sites)
c.dS.mutations <- calc.cumulative.muts(gene.dS.mutation.data, dS.normalization.const)

## Compare dN, dS, and  in the whole population!
c.total.dN.dS.plot <- plot.cumulative.muts(c.dN.mutations, my.color="purple", logscale=TRUE) %>%
add.cumulative.mut.layer(c.dS.mutations, my.color="green", logscale=TRUE)    

## add a layer for nonsense + indel + sv.
c.dN.dS.nonsense.indel.sv.plot <- c.total.dN.dS.plot %>%
    add.cumulative.mut.layer(c.sv.indel.nonsense.gene.mutations, "red", logscale=TRUE)

c.dN.dS.nonsense.indel.sv.plot
ggsave(c.dN.dS.nonsense.indel.sv.plot,filename="../results/figures/dN-dS-nonsenseindelsv.pdf")

## IMPORTANT TODO: look at Ben Good's code, and modify to calculate the number of missense,
## synonymous, and nonsense sites for normalization.


## As a comparison, look at distribution of mutations in intergenic regions.
## use this to study regulatory evolution, or perhaps relaxed selection?
intergenic.mutation.data <- filter(mutation.data, Gene=='intergenic')
intergenic.snp.data <- filter(intergenic.mutation.data,Annotation=='noncoding')
intergenic.sv.data <- filter(intergenic.mutation.data, Annotation=='sv')
intergenic.indel.data <- filter(intergenic.mutation.data, Annotation=='indel')

c.intergenic.mutations <- calc.cumulative.muts(intergenic.mutation.data,intergenic.length)
c.intergenic.snp.mutations <- calc.cumulative.muts(intergenic.snp.data,intergenic.length)
c.intergenic.sv.mutations <- calc.cumulative.muts(intergenic.sv.data,intergenic.length)
c.intergenic.indel.mutations <- calc.cumulative.muts(intergenic.indel.data,intergenic.length)

intergenic.plot <- plot.cumulative.muts(c.intergenic.snp.mutations,logscale=TRUE) %>%
    add.cumulative.mut.layer(c.intergenic.sv.mutations, "orange", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.intergenic.indel.mutations, "light blue", logscale=TRUE)

intergenic.plot

#############################################################################
## examining the rates that genes in the left-hand tail and right-hand tail
## of the positive selection distribution get hit by mutations.
## TODO: what is the relationship, if any, between genes with dN (or not)
## and those with nonsense.indels.sv's, (or not), and rates of accumulation?
## Is there a relationship between genes under positive or purifying selection?

### make plots for the top, median, and bottom genes by multinomial test.
## IMPORTANT NOTE: THE UNDERLYING SETS OF GENES NEED TO BE THE SAME.

#### Plots for multinom1to50 genes.
c.multinom1to50.dN <- filter(multinom.1to50.mutation.data, Annotation=='missense') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.dS <- filter(multinom.1to50.mutation.data, Annotation=='synonymous') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.nonsense <- filter(multinom.1to50.mutation.data, Annotation=='nonsense') %>%
    calc.cumulative.muts(multinom.1to50.length)
c.multinom1to50.nonsense.indel.sv <- filter(multinom.1to50.mutation.data,
                                           Annotation %in% c("nonsense","indel","sv")) %>%
    calc.cumulative.muts(multinom.1to50.length)

## Now make the plot.
log.multinom1to50.plot <- plot.cumulative.muts(c.multinom1to50.dN, my.color="purple") %>%
    add.cumulative.mut.layer(c.multinom1to50.dS, my.color="green") %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense, my.color="gray") %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense.indel.sv, my.color="black")

log.multinom1to50.plot
ggsave(multinom1to50.plot,filename="../results/figures/log-multinom-1to50.pdf")

multinom1to50.plot <- plot.cumulative.muts(c.multinom1to50.dN, my.color="purple",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.dS, my.color="green",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense, my.color="gray",logscale=FALSE) %>%
    add.cumulative.mut.layer(c.multinom1to50.nonsense.indel.sv, my.color="black",logscale=FALSE)

log.multinom1to50.plot
ggsave(multinom1to50.plot,filename="../results/figures/multinom-1to50.pdf")


##########################################################################
## can I "train" a model of relaxed selection on genes that we know are not under selection
## in the LTEE? That is, metabolic operons that are never used?

## Then perhaps I can use this model to dientangle positive selection from
## relaxed selection (if I am lucky.)


## I could use Brian Wade's genomes as a test set to validate
## predictions of genes under purifying selection (as a paper 3).

## based on redoing the binomial test on genomes, it seems the difference in result
## is NOT driven by target size. rather, it is whether or not all kinds of mutations
## are included, and not just restricting to point mutations.

##########################################################################
## look at accumulation of stars over time for genes in the different proteome
## sectors.
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN in each population.

## get proteome sector assignments from Hui et al. 2015 Supplementary Table 2.
## I saved a reduced version of the data.
proteome.assignments <- read.csv('../data/Hui-2015-proteome-section-assignments.csv',as.is=TRUE)
REL606.proteome.assignments <- inner_join(REL606.genes,proteome.assignments)

## add proteome assignment to gene mutation.data.
sector.mut.data <- inner_join(gene.mutation.data,REL606.proteome.assignments)

##six sectors:  "A" "S" "O" "U" "R" "C"
A.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='A')
S.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='S')
O.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='O')
U.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='U')
R.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='R')
C.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='C')

A.length <- sum(A.sector.mut.data$gene_length)
S.length <- sum(S.sector.mut.data$gene_length)
O.length <- sum(O.sector.mut.data$gene_length)
U.length <- sum(U.sector.mut.data$gene_length)
R.length <- sum(R.sector.mut.data$gene_length)
C.length <- sum(C.sector.mut.data$gene_length)

c.A.muts <- calc.cumulative.muts(A.sector.mut.data, A.length)
c.S.muts <- calc.cumulative.muts(S.sector.mut.data, S.length)
c.O.muts <- calc.cumulative.muts(O.sector.mut.data, O.length)
c.U.muts <- calc.cumulative.muts(U.sector.mut.data, U.length)
c.R.muts <- calc.cumulative.muts(R.sector.mut.data, R.length)
c.C.muts <- calc.cumulative.muts(C.sector.mut.data, C.length)

log.sector.plot <- plot.cumulative.muts(c.A.muts, my.color="black", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.mutations,my.color="grey", logscale=TRUE)
ggsave(log.sector.plot,filename='../results/figures/log-sector-plot.pdf')

sector.plot <- plot.cumulative.muts(c.A.muts, my.color="black", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.mutations,my.color="grey", logscale=FALSE)
ggsave(sector.plot,filename='../results/figures/sector-plot.pdf')

##########################################################################
## look at accumulation of stars over time for genes in different eigengenes
## inferred by Wytock and Motter (2018).
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum in each population.

## get eigengene sector assignments from Wytock and Motter (2018) Supplementary File 1.
## I saved a reduced version of the data.

eigengenes <- read.csv('../data/Wytock2018-eigengenes.csv',as.is=TRUE)
REL606.eigengenes <- inner_join(REL606.genes,eigengenes)

## add eigengene assignment to gene.mutation.data.
eigengene.mut.data <- inner_join(gene.mutation.data, REL606.eigengenes)

eigengene1.mut.data <- filter(eigengene.mut.data,Eigengene==1)
eigengene2.mut.data <- filter(eigengene.mut.data,Eigengene==2)
eigengene3.mut.data <- filter(eigengene.mut.data,Eigengene==3)
eigengene4.mut.data <- filter(eigengene.mut.data,Eigengene==4)
eigengene5.mut.data <- filter(eigengene.mut.data,Eigengene==5)
eigengene6.mut.data <- filter(eigengene.mut.data,Eigengene==6)
eigengene7.mut.data <- filter(eigengene.mut.data,Eigengene==7)
eigengene8.mut.data <- filter(eigengene.mut.data,Eigengene==8)
eigengene9.mut.data <- filter(eigengene.mut.data,Eigengene==9)

eigen1.length <- sum(eigengene1.mut.data$gene_length)
eigen2.length <- sum(eigengene2.mut.data$gene_length)
eigen3.length <- sum(eigengene3.mut.data$gene_length)
eigen4.length <- sum(eigengene4.mut.data$gene_length)
eigen5.length <- sum(eigengene5.mut.data$gene_length)
eigen6.length <- sum(eigengene6.mut.data$gene_length)
eigen7.length <- sum(eigengene7.mut.data$gene_length)
eigen8.length <- sum(eigengene8.mut.data$gene_length)
eigen9.length <- sum(eigengene9.mut.data$gene_length)

c.eigen1.muts <- calc.cumulative.muts(eigengene1.mut.data, eigen1.length)
c.eigen2.muts <- calc.cumulative.muts(eigengene2.mut.data, eigen2.length)
c.eigen3.muts <- calc.cumulative.muts(eigengene3.mut.data, eigen3.length)
c.eigen4.muts <- calc.cumulative.muts(eigengene4.mut.data, eigen4.length)
c.eigen5.muts <- calc.cumulative.muts(eigengene5.mut.data, eigen5.length)
c.eigen6.muts <- calc.cumulative.muts(eigengene6.mut.data, eigen6.length)
c.eigen7.muts <- calc.cumulative.muts(eigengene7.mut.data, eigen7.length)
c.eigen8.muts <- calc.cumulative.muts(eigengene8.mut.data, eigen8.length)
c.eigen9.muts <- calc.cumulative.muts(eigengene9.mut.data, eigen9.length)

log.eigen.plot <- plot.cumulative.muts(c.eigen1.muts, my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.mutations,my.color="grey", logscale=TRUE)
ggsave(log.eigen.plot,filename='../results/figures/log-eigen-plot.pdf')

eigen.plot <- plot.cumulative.muts(c.eigen1.muts, my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.mutations,my.color="grey", logscale=FALSE)
ggsave(eigen.plot,filename='../results/figures/eigen-plot.pdf')

##########################################################################
## Bootstrap a distribution for the distribution of the accumulation of stars
## over time for random subsets of genes.
## Say, 1000 or 10000 random subsets of genes.
## This idea could be used to bootstrap p-values for statistics, by
## comparing subsets of mutations in genes of interest to this random null model.
## Distributions greater or lower than all 1000 curves is significantly greater or lesser
## at a p = 0.001/2 (I think? check this calculation more rigorously.)

plot.random.subsets <- function(data, subset.size=300, N=1000,log=TRUE) {

  ## set up an empty plot then add random trajectories, one by one.
  my.plot <- ggplot(data) +
  theme_classic() +
  facet_wrap(.~Population,scales='fixed') +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')

  for (i in 1:N) {
    rando.genes <- sample(unique(data$Gene),subset.size)
    mut.subset <- filter(data,Gene %in% rando.genes)
    subset.length <- sum(mut.subset$gene_length)
    c.mut.subset <- calc.cumulative.muts(mut.subset,subset.length)
    
    if (log) {
      my.plot <- my.plot +
      geom_point(data=c.mut.subset,aes(x=Generation,y=log(normalized.cs)), color='gray',size=0.2,alpha = 0.1)
    } else {
      my.plot <- my.plot +
      geom_point(data=c.mut.subset,aes(x=Generation,y=normalized.cs), color='gray',size=0.2,alpha = 0.1)
    }
  }
  return(my.plot)
}

log.rando.plot <- plot.random.subsets(gene.mutation.data, log=TRUE)
ggsave(log.rando.plot,filename='../results/figures/log-rando-plot.png')
rando.plot <- plot.random.subsets(gene.mutation.data,log=FALSE)
ggsave(rando.plot,filename='../results/figures/rando-plot.png')

#######################################################################################
## TODO: infer cohorts using Haixu Tang's new algorithm.
## Then ask whether cohorts show functional enrichment using
## STRING annotation.
