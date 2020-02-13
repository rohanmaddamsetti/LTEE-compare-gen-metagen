## metagenomics-analysis.R by Rohan Maddamsetti.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## TODO: measure the concentration/density of mutation as an empirical CDF.

## TODO: MAKE SURE MASKED REGIONS ARE TREATED PROPERLY!

## TODO: It could be quite interesting to examine genes in prophage that
## Jeff Barrick has annotated as misc_regions.

library(tidyverse)

##########################################################################
## FUNCTIONS
##########################################################################
## look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.
## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length).
calc.cumulative.muts <- function(data, normalization.constant=NA) {

    cumsum.per.pop.helper.func <- function(df) {
        df %>%
            arrange(t0) %>%
            group_by(Population,Generation) %>%
            summarize(count=n()) %>%
            mutate(cs=cumsum(count)) %>%
            ungroup() 
    }
    
    ## if normalization.constant is not provided, then
    ## calculate based on gene length by default.
    if (is.na(normalization.constant)) {
        my.genes <- data %>% select(Gene,gene_length) %>% distinct()
        normalization.constant <- sum(my.genes$gene_length)
    }
    
    c.dat <- data %>%
        split(.$Population) %>%
        map_dfr(.f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    return(c.dat)
}

## calculate the tail probabilities of the true cumulative mutation trajectory
## of a given vector of genes (a 'module'), based on resampling
## random sets of genes. Returns both upper tail of null distribution,
## or P(random trajectory >= the actual trajectory).
## Output: a dataframe with three columns: Population, count, p.val
calculate.trajectory.tail.probs <- function(data, gene.vec, N=10000, normalization.constant=NA) {

    ## resamples have the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)
    
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- sample(unique(data$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset, normalization.constant) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    gene.vec.data <- data %>% filter(Gene %in% gene.vec)
    data.trajectory <- calc.cumulative.muts(gene.vec.data,normalization.constant)
    data.trajectory.summary <- data.trajectory %>%
        group_by(Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup()
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 

    trajectory.filter.helper <- function(pop.trajectories) {
        pop <- unique(pop.trajectories$Population)
        data.traj <- filter(data.trajectory.summary,Population == pop)
        final.data.norm.cs <- unique(data.traj$final.norm.cs)
        tail.trajectories <- filter(pop.trajectories, final.norm.cs >= final.data.norm.cs)
        return(tail.trajectories)
    }
    
    ## split by Population, then filter for bootstraps > data trajectory.
    uppertail.probs <- trajectory.summary %>%
        split(.$Population) %>%
        map_dfr(.f=trajectory.filter.helper) %>%
        group_by(Population) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
        
    return(uppertail.probs)
}

## Calculate the derivative of the cumulative accumulation of mutation occurrence.
## This is simply the rate of mutation occurrence in a class of genes.
calc.slope.of.cumulative.muts <- function(c.muts) {

    calc.slope.per.pop.helper.func <- function(df) {
        df %>%
            group_by(Population) %>%
            mutate(D.cs = cs - lag(cs)) %>%
            mutate(D.normalized.cs = normalized.cs - lag(normalized.cs))
    }
    
    D.of.c.muts <- c.muts %>%
        split(.$Population) %>%
        map_dfr(.f=calc.slope.per.pop.helper.func) %>%
        ## remove any NA values.
        na.omit()
    return(D.of.c.muts)
}


## This plot visualizes a two-tailed test (alpha = 0.05)
## against a bootstrapped null distribution.
## This is what we want to use for publication.
plot.base.layer <- function(data, subset.size=300, N=1000, alpha = 0.05, logscale=FALSE, normalization.constant=NA) {

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- sample(unique(data$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset, normalization.constant) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    ## filter out the top alpha/2 and bottom alpha/2 trajectories from each population,
    ## for a two-sided test. default is alpha == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(alpha/2) %>%
        select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(-alpha/2) %>%
        select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        select(-in.top,-in.bottom)

    if (logscale) {
        p <- ggplot(filtered.trajectories,aes(x=Generation,y=log10(normalized.cs))) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(filtered.trajectories,aes(x=Generation,y=normalized.cs)) +
            ylab('Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color='gray') +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14))
    return(p)                
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=FALSE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), size=0.2, color=my.color)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, aes(x=Generation,y=normalized.cs), size=0.2, color=my.color)
        }
    return(p)
}

## calculate cumulative numbers of mutations in each category.
## for vanilla plotting, without null distributions, as plotted by
## plot.base.layer.
plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log10(normalized.cs))) +
            ##ylim(-7,-2) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs)) +
            ##ylim(0,0.003) +
            ylab('Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)
    return(p)
}

## calculate derivative of cumulative numbers of mutations in each category.
plot.slope.of.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log10(D.normalized.cs))) +
            ylab('log[Slope of Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=D.normalized.cs)) +
            ylab('Slope of Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        geom_smooth(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free') +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)
    return(p)
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.slope.of.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=TRUE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log10(D.normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, aes(x=Generation,y=log10(D.normalized.cs)), color=my.color, size=0.2) +
            geom_smooth(size=0.2, color=my.color)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=D.normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, aes(x=Generation,y=D.normalized.cs), color=my.color, size=0.2) +
                geom_smooth(size=0.2, color=my.color)
        }
    return(p)
}

########################################################################

ks.analysis <- function(the.data, REL606.genes) {
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
    the.results.to.plot$index <- seq_len(nrow(the.results.to.plot))
  
    plot <- ggplot(the.results.to.plot, aes(x=index)) +
        geom_line(aes(y=empirical), colour="red") + 
        geom_line(aes(y=null), linetype=2) + 
        scale_x_continuous('Genes ranked by length',limits=c(0,4400)) +
        scale_y_continuous('Cumulative proportion of mutations',limits=c(0,1)) +
        theme_classic() +
        theme(axis.title=element_text(size=18),axis.text=element_text(size=12))
    
    return(plot)
}

################################################################################
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

    ## CRITICAL STEP: replace NAs with zeros.
    ## We need to keep track of genes that haven't been hit by any mutations
    density.df[is.na(density.df)] <- 0
    density.df <- tbl_df(density.df)

    
    return(density.df)
}

## Examine the gene mutation density by population.
pop.calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    density.df <- gene.mutation.data %>%
        filter(Annotation %in% mut_type_vec) %>%
        filter(Gene!= "intergenic") %>%
        mutate(Gene=as.factor(Gene)) %>%
        group_by(Population, Gene, gene_length) %>%
        summarize(mut.count=n()) %>%
        ungroup() %>%
        mutate(density=mut.count/gene_length) %>%
        arrange(desc(density))

    ## CRITICAL STEP: replace NAs with zeros.
    ## We need to keep track of genes that haven't been hit by any mutations
    density.df[is.na(density.df)] <- 0
    density.df <- tbl_df(density.df)
    
    return(density.df)
}

##########################################################################
## DATA ANALYSIS
##########################################################################

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

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops)))


gene.mutation.data <- inner_join(mutation.data,REL606.genes)

#############!!!!!!!!!!!!!!!!!!!!!!!!!
### IMPORTANT: THIS IS AN OBVIOUS SOURCE OF BUGS--
## double-check whether intergenic mutations
## are included or not-- as appropriate-- in all references to c.mutations
## here and in the aerobic/anaerobic code.

## for now, ONLY use mutation.data when examining mutation bias over the genome.
## when examining evolution in different gene sets, use gene.mutation.data.


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


## let's combine IS (structural mutations), indels, and nonsense mutations,
## as a proxy for purifying selection.
gene.nonsense.sv.indels.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation %in% c('indel', 'sv', 'nonsense'))

## Now filter for nonsense, sv, and indels in genes (of course nonsense are always in genes.)
sv.indel.nonsense.gene.mutation.data <- gene.mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))

## filter for point mutations (missense, synonymous, nonsense).
gene.point.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation %in% c('missense', 'synonymous', 'nonsense'))

## look at distribution over the genome for different classes of mutations.

gene.only.mutation.data <- gene.mutation.data %>%
    filter(Annotation !='noncoding')

gene.dS.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='synonymous')

gene.dN.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='missense')

gene.nonsense.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='nonsense')

gene.noncoding.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='noncoding')

gene.except.dS.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation!='synonymous')

gene.sv.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation=='sv')

gene.indel.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation=='indel')
##########################################################################################
## make base plot of null distributions by subsampling/bootstrapping.
## I will overlay data on top of these plots.

## Base plots here of null distributions: add the data lines on top to compare.
log.all.rando.plot <- plot.base.layer(gene.mutation.data,logscale=TRUE)

log.sv.indel.nonsen.rando.plot <- plot.base.layer(sv.indel.nonsense.gene.mutation.data, logscale=TRUE)

all.rando.plot <- plot.base.layer(gene.mutation.data, logscale=FALSE)
sv.indel.nonsen.rando.plot <- plot.base.layer(sv.indel.nonsense.gene.mutation.data, logscale=FALSE)

##########################################################################################
## investigate dS (ano other classes of mutations) across the genome
## in the metagenomics data.
## revamp code from my 2015 Mol. Biol. Evol. paper.
## in short, cannot reject null model that dS is uniform over the genome.
## and dN fits the null extremely well-- even better than dS!
## Probably because there are 3 times as many dN as dS throughout the
## experiment.

## examine all mutations over genes in the genome.
cumsum.all.over.metagenome <- ks.analysis(gene.only.mutation.data,REL606.genes)
all.KS.plot <- make.KS.Figure(cumsum.all.over.metagenome)

## examine dS over the genome.
cumsum.dS.over.metagenome <- ks.analysis(gene.dS.mutation.data,REL606.genes)
dS.KS.plot <- make.KS.Figure(cumsum.dS.over.metagenome)

## examing dN over the genome.

## dN fits the null even better than dS! But note that
## there are 18493 dN in the data, and 6792 dS in the data.
## so the better fit is probably best explained by the larger sample size.
cumsum.dN.over.metagenome <- ks.analysis(gene.dN.mutation.data,REL606.genes)
make.KS.Figure(cumsum.dN.over.metagenome)

## let's look at nonsense mutations.
## nonsense mutations don't fit the null expectation.
## opposite trend of indels or IS elements, though!
## not sure why.
cumsum.nonsense.over.metagenome <- ks.analysis(gene.nonsense.mutation.data,REL606.genes)
make.KS.Figure(cumsum.nonsense.over.metagenome)

## examine non-coding mutations.
## very clear that the normalization by gene_length is not appropriate.
cumsum.noncoding.over.metagenome <- ks.analysis(gene.noncoding.mutation.data,REL606.genes)
make.KS.Figure(cumsum.noncoding.over.metagenome)

## now let's look at all mutations except for dS.
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

## Nonsense mutation distribution is predicted by neither the distribution
## of synonymous nor missense mutations! Nice result.
ks.test(cumsum.dS.over.metagenome$empirical,
        cumsum.nonsense.over.metagenome$empirical,
        simulate.p.value=TRUE)

ks.test(cumsum.dN.over.metagenome$empirical,
        cumsum.nonsense.over.metagenome$empirical,
        simulate.p.value=TRUE)


## VERY IMPORTANT TODO: EXCLUDE SV AND INDELS THAT ARE ANNOTATED BY A GENE,
## BUT ACTUALLY OCCUR IN PROMOTER REGIONS!

## examine structural mutations (IS elements) affecting genes.
## RESULT: longer genes are depleted in IS insertions!!
cumsum.sv.over.metagenome <- ks.analysis(gene.sv.mutation.data,REL606.genes)
make.KS.Figure(cumsum.sv.over.metagenome)

## examine indels.
## RESULT: longer genes are depleted in indels!
cumsum.indel.over.metagenome <- ks.analysis(gene.indel.mutation.data,REL606.genes)
make.KS.Figure(cumsum.indel.over.metagenome)

## Therefore, structural mutations and indels are probably what are driving the
## differences that I saw in aerobic and anerobic mutations before.

cumsum.nonsense.sv.indels.over.metagenome <- ks.analysis(gene.nonsense.sv.indels.mutation.data, REL606.genes)
make.KS.Figure(cumsum.nonsense.sv.indels.over.metagenome)

## Figure 1 will be a conceptual figure of the different evolutionary models
## for different sets of genes.
## (purifying selection, coupon collecting, neutral or mutation accumulation, etc.)

## Figure 2 TODO: 
## combine these results to show distribution of different classes of mutations
## over the genome.

##########################################################################################
## MUTATION BIAS ANALYSIS.

## TODO: plot distribution of structural variants and indels over the genome for each
## population.

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

## plot the distribution of indels and SV over the genome, separately for the
## different LTEE populations, and altogether.


##########################################################################################
## cross-check the KS-test results with 50K genomes.
## These were downloaded from Jeff Barrick's shiny app:
## https://barricklab.org/shiny/LTEE-Ecoli.

gene.50K.mutations <- read.csv("../data/Gen50000_allMutations.csv",
                              header=TRUE,
                              as.is=TRUE) %>%
    filter(clone=='A') %>%
    rename(Gene=gene_name) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

dS.50K <- filter(gene.50K.mutations, snp_type=='synonymous')
no.dS.50K <- filter(gene.50K.mutations, snp_type!='synonymous')
dN.50K <- filter(gene.50K.mutations, snp_type!='nonsynonymous')
IS.50K <- filter(gene.50K.mutations, mutation_category=="mobile_element_insertion")
indels.50K <- filter(gene.50K.mutations, mutation_category=="small_indel")
nonsense.50K <- filter(gene.50K.mutations, mutation_category=="snp_nonsense")
largedeletions.50K <- filter(gene.50K.mutations, mutation_category=="large_deletion")

## cross-check with dS in 50K genomes.
cumsum.dS.genome.50K <- ks.analysis(dS.50K, REL606.genes)
make.KS.Figure(cumsum.dS.genome.50K)

## cross-check with everything except dS in 50K genomes.
cumsum.no.dS.genome.50K <- ks.analysis(no.dS.50K, REL606.genes)
make.KS.Figure(cumsum.no.dS.genome.50K)

## cross-check with dN in 50K genomes.
cumsum.dN.genome.50K <- ks.analysis(dN.50K, REL606.genes)
make.KS.Figure(cumsum.dN.genome.50K)

## cross-check IS elements in the 50K genomes.
cumsum.IS.genome.50K <- ks.analysis(IS.50K,REL606.genes)
make.KS.Figure(cumsum.IS.genome.50K)

## cross-check indels in the 50K genomes.
cumsum.indels.genome.50K <- ks.analysis(indels.50K,REL606.genes)
make.KS.Figure(cumsum.indels.genome.50K)

## cross-check nonsense SNPs in the 50K genomes.
## OPPOSITE trend!!! Not sure how to interpret...
cumsum.nonsense.genome.50K <- ks.analysis(nonsense.50K,REL606.genes)
make.KS.Figure(cumsum.nonsense.genome.50K)

## large deletions in 50K genomes.
cumsum.largedeletions.genome.50K <- ks.analysis(largedeletions.50K,REL606.genes)
make.KS.Figure(cumsum.largedeletions.genome.50K)

#################################################################################
### IMPORTANT BUG! GENES THAT SHOULD HAVE BEEN MASKED, WITH NO
### MUTATION CALLS IN THE GOOD ET AL. DATA, REAPPEAR IN REL606 GENES!!!
## PROPHAGE STARTING AT ECB_00814 is a good place to start debugging.
################################################################################

## Let's look at mutation densities within populations (not summed over populations.)

## add an index column to the dataframe.
add_index_helper <- function(df) {
    df %>% mutate(index = seq_len(nrow(.)))
}

pop.all.mutation.density <- gene.mutation.data %>%
    pop.calc.gene.mutation.density(mut_type_vec=c("missense", "sv", "synonymous","noncoding", "indel", "nonsense")) %>% arrange(Population,desc(density)) %>%
    split(.$Population) %>%
    map_dfr(add_index_helper)

pop.all.mutation.density.plot <- ggplot(pop.all.mutation.density,
                                        aes(x=index,y=density)) +
    geom_point(size=0.3) + theme_classic() + facet_wrap(.~Population,scale='fixed',nrow=4)

pop.all.mutation.density.plot

## look at dN now.
pop.dN.mutation.density <- gene.mutation.data %>%
    pop.calc.gene.mutation.density(mut_type_vec=c("missense")) %>% arrange(Population,desc(density)) %>%
    split(.$Population) %>%
    map_dfr(add_index_helper)

pop.dN.mutation.density.plot <- ggplot(pop.dN.mutation.density,
                                        aes(x=index,y=density)) +
    geom_point(size=0.3) + theme_classic() + facet_wrap(.~Population,scale='fixed',nrow=4)

pop.dN.mutation.density.plot


################################################################################
### DENSITIES SUMMED OVER ALL POPULATIONS.

## IMPORTANT: these densities are summed over ALL LTEE populations,
## so they don't correspond to single LTEE populations.

## just rank gene.mutation.density to see what it looks like!
## this is already a good indication of selection,
## and this is pretty much what Ben Good does in his Table S3.

## IMPORTANT: these densities are summed over ALL LTEE populations,
## so they don't correspond to single LTEE populations.

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

nonsense.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data, c("nonsense")) %>%
    rename(nonsense.mut.count = mut.count) %>%
    rename(nonsense.mut.density = density)

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
    full_join(nonsense.mutation.density) %>%
    full_join(sv.mutation.density) %>%
    full_join(indel.mutation.density) %>%
    full_join(nonsense.indel.sv.density) %>%
    full_join(all.except.dS.density)

#### CRITICAL STEP: replace NAs with zeros.
#### We need to keep track of genes that haven't been hit by any mutations
#### in a given mutation class (sv, indels, dN, etc.)
gene.mutation.densities[is.na(gene.mutation.densities)] <- 0
gene.mutation.densities <- tbl_df(gene.mutation.densities)

##################################
## break REL606 genes into deciles, based on mutation density.

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


## seems like a lot of the all.density.data results are driven by repeated IS
## insertions.

## TODO: plot IS insertions over the genome, and indels over the genome, when
## examining mutation bias.

## let's look at nonsense, sv, and indels, and see what the distribution looks like.
nonindelsv.density.data <- arrange(gene.mutation.densities,desc(nonsense.indel.sv.mut.density)) %>%     ## for plotting convenience, add an index to the data frame.
    mutate(index = seq_len(nrow(.)))

nonindelsv.density.plot <- ggplot(nonindelsv.density.data, aes(x=index,y=nonsense.indel.sv.mut.density)) +
    geom_point() + theme_classic()

nonindelsv.density.plot

## 2,735 genes have no nonsense, indels, or structural variants throughout the LTEE.
zero.nonindelsv.density.data <- nonindelsv.density.data %>%
    filter(nonsense.indel.sv.mut.density==0)

top.nonindelsv.density.data <- nonindelsv.density.data %>%
    filter(nonsense.indel.sv.mut.density>0.01) %>%
    select(Gene,locus_tag,blattner,gene_length,product,nonsense.indel.sv.mut.density,nonsense.indel.sv.mut.count)



## let's just look at dN and see what the distribution looks like.
dN.density.data <- arrange(gene.mutation.densities,desc(dN.mut.density)) %>%
    ## for plotting convenience, add an index to the data frame.
    mutate(index = seq_len(nrow(.))) %>%
    mutate(top=index<=420)

dN.density.plot <- ggplot(dN.density.data, aes(x=index,y=dN.mut.density,color=top)) +
    geom_point() + theme_classic()

dN.density.plot

## let's quickly plot cumulative numbers of mutations in the top genes
## (damn the circularity, full speed ahead!)

top.dN.gene.mutation.data <- filter(gene.dN.mutation.data,Gene %in% top.by.dN.density.genes$Gene)

c.top.dN.gene.mutations <- calc.cumulative.muts(top.dN.gene.mutation.data)

log.c.top.dN.gene.plot <- log.small.dN.rando.plot %>%
    add.cumulative.mut.layer(c.top.dN.gene.mutations, my.color="black",logscale=TRUE)

c.top.dN.gene.plot <- small.dN.rando.plot %>%
    add.cumulative.mut.layer(c.top.dN.gene.mutations, my.color="black",logscale=FALSE)

## what are these top genes?
top.by.dN.density.genes

## let's filter further on dN.mut.density.
top_dN_checkme <- top.by.dN.density.genes %>% select(Gene,locus_tag,blattner, gene_length,
                                                 product, dN.mut.count,dN.mut.density,
                                                 dN.decile) %>%
    filter(dN.mut.density>0.01) %>%
    arrange(desc(dN.mut.density))
## Curious list of genes... let's check STRING.
write.csv(top_dN_checkme,file="../results/dN_checkme.csv")

## Probably better to import gene lists from Tenaillon and Good papers to avoid the
## obvious circularity. But EVEN better to compare these lists. When are they the same;
## when are they different? Where do the genes from Tenaillon and Good papers
## fall in the dN density distribution?

## interesting: ribosomal genes are in both the top as well as the bottom
## (purifying selection results).

top.dN.ribosome.genes <- c('yceD', 'yeiP', 'rplL', 'rpsD', 'rpsE', 'infB', 'rplF', 'rplS',
                           'rplP', 'rplY', 'ykgM')
top.dN.ribosome.mutations <- filter(gene.dN.mutation.data,Gene %in% top.dN.ribosome.genes)
c.top.dN.ribosome.mutations <- calc.cumulative.muts(top.dN.ribosome.mutations)

log.c.top.dN.ribosome.plot <- log.small.dN.rando.plot %>%
    add.cumulative.mut.layer(c.top.dN.ribosome.mutations, my.color="black",logscale=TRUE)

c.top.dN.ribosome.plot <- small.dN.rando.plot %>%
    add.cumulative.mut.layer(c.top.dN.ribosome.mutations, my.color="black",logscale=FALSE)

## quick random check: look for nucleotide-level parallelism in the Good dataset.
nuc.parallel.data <- gene.mutation.data %>% group_by(Gene, Position, Annotation, gene_length, product, blattner, locus_tag) %>% summarize(count=n()) %>%
    filter(count>=3) %>% arrange(desc(count))

write.csv(nuc.parallel.data,'../results/parallel-nuc-summary.csv')

#######################################################################################
## INITIAL PURIFYING SELECTION RESULTS.
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
##target.size.numbers <- read.csv('../results/target_size.csv',header=TRUE,as.is=TRUE)

##total.length <- filter(target.size.numbers,set=='genome')$total_gene_length
##total.synon.sites <- filter(target.size.numbers,set=='genome')$synon_sites
##total.nonsynon.sites <- filter(target.size.numbers,set=='genome')$non_synon_sites

######################NOTE: DOUBLE CHECK CONSISTENT WITH measureTargetSize.py. output!!!!!
## IMPORTANT TODO: THESE NUMBERS ARE ALL OFF NOW-- HAVE TO TAKE MASKING INTO ACCOUNT!!!
####### Constants. REVAMP THIS CODE!
## from measureIntergenicTargetSize.py:
## Length of intergenic regions: 487863
##intergenic.length <- 487863


########################
## look at the accumulation of all sv, indels, and nonsense mutations
## in each population.

c.sv.indel.nonsense.gene.mutations <- calc.cumulative.muts(sv.indel.nonsense.gene.mutation.data)
c.sv.indel.nonsense.gene.plot <- sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.sv.indel.nonsense.gene.mutations,my.color='black')


########################
## Examine ALL mutations, including intergenic regions.
#############!!!!!!!!!!!!!!!!!!!!!!!!!
### IMPORTANT: THIS IS AN OBVIOUS SOURCE OF BUGS--
## double-check whether intergenic mutations
## are included or not-- as appropriate-- in all references to c.mutations
## here and in the aerobic/anaerobic code.
c.mutations <- calc.cumulative.muts(gene.mutation.data)

##########################
## let's calculate the rate of mutation occurrence
## (this is the derivative of the cumulative number of mutations over time).

D.of.c.mutations <- calc.slope.of.cumulative.muts(c.mutations)
D.of.c.mutations.plot <- plot.slope.of.cumulative.muts(D.of.c.mutations, logscale=TRUE)

D.of.c.mutations.plot

c.mutation.plot <- plot.cumulative.muts(c.mutations, logscale=TRUE)
c.mutation.plot2 <- plot.cumulative.muts(c.mutations, logscale=FALSE)

## Examine dN mutations.
##c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
##dN.normalization.const <- total.length * total.nonsynon.sites/(total.synon.sites+total.nonsynon.sites)
c.dN.mutations <- calc.cumulative.muts(gene.dN.mutation.data)

## Examine dS mutations.
##c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)
##dS.normalization.const <- total.length * total.synon.sites/(total.synon.sites+total.nonsynon.sites)
c.dS.mutations <- calc.cumulative.muts(gene.dS.mutation.data)

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

c.A.muts <- calc.cumulative.muts(A.sector.mut.data)
c.S.muts <- calc.cumulative.muts(S.sector.mut.data)
c.O.muts <- calc.cumulative.muts(O.sector.mut.data)
c.U.muts <- calc.cumulative.muts(U.sector.mut.data)
c.R.muts <- calc.cumulative.muts(R.sector.mut.data)
c.C.muts <- calc.cumulative.muts(C.sector.mut.data)

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

## just plot sv, indels, and nonsense mutations.

A.sv.indel.nonsen.muts <- filter(A.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
S.sv.indel.nonsen.muts <- filter(S.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
O.sv.indel.nonsen.muts <- filter(O.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
U.sv.indel.nonsen.muts <- filter(U.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
R.sv.indel.nonsen.muts <- filter(R.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
C.sv.indel.nonsen.muts <- filter(C.sector.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
    
c.A.sv.indel.nonsen.muts <- calc.cumulative.muts(A.sv.indel.nonsen.muts)
c.S.sv.indel.nonsen.muts <- calc.cumulative.muts(S.sv.indel.nonsen.muts)
c.O.sv.indel.nonsen.muts <- calc.cumulative.muts(O.sv.indel.nonsen.muts)
c.U.sv.indel.nonsen.muts <- calc.cumulative.muts(U.sv.indel.nonsen.muts)
c.R.sv.indel.nonsen.muts <- calc.cumulative.muts(R.sv.indel.nonsen.muts)
c.C.sv.indel.nonsen.muts <- calc.cumulative.muts(C.sv.indel.nonsen.muts)


## plot cumulative plots for proteome sectors.
log.sector.testplot <- log.all.rando.plot %>%
    add.cumulative.mut.layer(c.A.muts, my.color="black", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange", logscale=TRUE)
ggsave(log.sector.testplot,filename='../results/figures/log-sector-testplot.png')

sector.testplot <- all.rando.plot %>%
    add.cumulative.mut.layer(c.A.muts, my.color="black", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange", logscale=FALSE)
ggsave(sector.testplot,filename='../results/figures/sector-testplot.png')


## now plot sv, indel, nonsense mutations on top.
log.sector.sv.indel.nonsen.testplot <- log.sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.A.sv.indel.nonsen.muts, my.color="black", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.S.sv.indel.nonsen.muts,my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.O.sv.indel.nonsen.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.U.sv.indel.nonsen.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.R.sv.indel.nonsen.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.C.sv.indel.nonsen.muts,my.color="orange", logscale=TRUE)
ggsave(log.sector.sv.indel.nonsen.testplot,filename='../results/figures/log-sector-sv-indel-nonsen-testplot.png')

sector.sv.indel.nonsen.testplot <- sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.A.sv.indel.nonsen.muts, my.color="black", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.S.sv.indel.nonsen.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.O.sv.indel.nonsen.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.U.sv.indel.nonsen.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.R.sv.indel.nonsen.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.C.sv.indel.nonsen.muts,my.color="orange", logscale=FALSE)
ggsave(sector.sv.indel.nonsen.testplot,filename='../results/figures/sector-sv-indel-nonsen-testplot.png')

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

c.eigen1.muts <- calc.cumulative.muts(eigengene1.mut.data)
c.eigen2.muts <- calc.cumulative.muts(eigengene2.mut.data)
c.eigen3.muts <- calc.cumulative.muts(eigengene3.mut.data)
c.eigen4.muts <- calc.cumulative.muts(eigengene4.mut.data)
c.eigen5.muts <- calc.cumulative.muts(eigengene5.mut.data)
c.eigen6.muts <- calc.cumulative.muts(eigengene6.mut.data)
c.eigen7.muts <- calc.cumulative.muts(eigengene7.mut.data)
c.eigen8.muts <- calc.cumulative.muts(eigengene8.mut.data)
c.eigen9.muts <- calc.cumulative.muts(eigengene9.mut.data)

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

## just plot sv, indels, and nonsense mutations.
sv.indel.nonsen.eigengene1.muts <- filter(eigengene1.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene2.muts <- filter(eigengene2.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene3.muts <- filter(eigengene3.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene4.muts <- filter(eigengene4.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene5.muts <- filter(eigengene5.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene6.muts <- filter(eigengene6.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene7.muts <- filter(eigengene7.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene8.muts <- filter(eigengene8.mut.data,Annotation %in% c("sv", "indel", "nonsense"))
sv.indel.nonsen.eigengene9.muts <- filter(eigengene9.mut.data,Annotation %in% c("sv", "indel", "nonsense"))

c.sv.indel.nonsen.eigen1.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene1.muts)
c.sv.indel.nonsen.eigen2.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene2.muts)
c.sv.indel.nonsen.eigen3.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene3.muts)
c.sv.indel.nonsen.eigen4.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene4.muts)
c.sv.indel.nonsen.eigen5.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene5.muts)
c.sv.indel.nonsen.eigen6.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene6.muts)
c.sv.indel.nonsen.eigen7.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene7.muts)
c.sv.indel.nonsen.eigen8.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene8.muts)
c.sv.indel.nonsen.eigen9.muts <- calc.cumulative.muts(sv.indel.nonsen.eigengene9.muts)



log.eigen.testplot <- log.all.rando.plot %>%
    add.cumulative.mut.layer(c.eigen1.muts,my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black", logscale=TRUE)
ggsave(log.eigen.testplot,filename='../results/figures/log-eigen-testplot.png')

eigen.testplot <- all.rando.plot %>%
    add.cumulative.mut.layer(c.eigen1.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black", logscale=FALSE)
ggsave(eigen.testplot,filename='../results/figures/eigen-testplot.png')


## plot eigen sv, indels, nonsense.
log.eigen.sv.indel.nonsen.testplot <- log.sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen1.muts,my.color="red", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen2.muts,my.color="orange", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen3.muts,my.color="yellow", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen4.muts,my.color="green", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen5.muts,my.color="cyan", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen6.muts,my.color="blue", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen7.muts,my.color="violet", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen8.muts,my.color="pink", logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen9.muts,my.color="black", logscale=TRUE)
ggsave(log.eigen.sv.indel.nonsen.testplot,filename='../results/figures/log-eigen-sv-indel-nonsen-testplot.png')

eigen.sv.indel.nonsen.testplot <- sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen1.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen2.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen3.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen4.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen5.muts,my.color="cyan", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen6.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen7.muts,my.color="violet", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen8.muts,my.color="pink", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen9.muts,my.color="black", logscale=FALSE)
ggsave(eigen.sv.indel.nonsen.testplot,filename='../results/figures/eigen-sv-indel-nonsen-testplot.png')

##########################################################################
## NOTES:
##########################################################################
## can I "train" a model of relaxed selection on genes that we know are not under selection
## in the LTEE? That is, metabolic operons that are never used?

## Then perhaps I can use this model to disentangle positive selection from
## relaxed selection (if I am lucky.)

## I think so! just bootstrap/resample trajectories from a training set of
## metabolic operons whose function decay in the LTEE from previous work.
## This could be a good baseline for comparing other sets of genes.
## (of course, result depends on how well this set is chosen.
## bootstrapping random sets from the whole genome may be a safer comparison).


## I could use Brian Wade's genomes as a test set to validate
## predictions of genes under purifying selection (as a paper 3).

## based on redoing the binomial test on genomes, it seems the difference in result
## is NOT driven by target size. rather, it is whether or not all kinds of mutations
## are included, and not just restricting to point mutations.

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

