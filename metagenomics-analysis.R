## metagenomics-analysis.R by Rohan Maddamsetti and Nkrumah Grant 

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## IMPORTANT TODO:
## The problem with setting up a probabilistic model for mutation rate variation is
## that causal factors evolve during the evolution experiment! Modeling these
## changes adds a lot of complexity. For this reason, rather than directly
## modeling mutation occurrence over the genome as a periodic function of
## distance from the replication origin, as as a function of leading/lagging strand,
## we used importance sampling to sample randomized modules of genes based on the empirical
## density of mutations per gene. We could also use this technique to sample genes
## using a sampling function that takes gene length, leading/lagging strand and a periodic
## function of distance from the oriC replication origin.

## If I need an independent dataset to set phase parameters-- then use E. coli MA
## experiment data to calibrate. Use data in Foster et al. (2013) in G3,
## or use any more recent datasets from that group or others.

##to calculate an intensity value
##over the genome (regression approach to disentangle, citing Bailey & Bataillon paper.

## NOTE: Perhaps should cite my STLE paper-- recombination events tend to happen flanking
## the replication origin. Is that result connected to the wave-like mutation bias
## pattern seen here and in Patricia Foster's and Vaughn Cooper's evolution experiments?
## Also check out Hi-C papers and others reporting 3D chromosome in E. coli.

library(tidyverse)
library(cowplot)

##########################################################################
## FUNCTIONS FOR DATA ANALYSIS REPORTED IN MAIN TEXT
##########################################################################
## calculate the probability that a locus (or set of loci) is not hit by mutations,
## assuming uniform mutation rate.
## l is locus length, and n is the number of mutations in the dataset.
probability.that.not.hit <- function(l,n) {
    GENOME.LENGTH <- 4629812
    p <- (1 - (l/GENOME.LENGTH))^n
    return(p)
}

#' function to rotate REL606 genome coordinates, setting oriC/terB at the center of plots
#' that examine mutation bias over the chromosome.
rotate.REL606.chr <- function(my.position, c) {
    ## we want to change coordinates so that c is the new origin.
    GENOME.LENGTH <- 4629812
    midpoint <- GENOME.LENGTH/2
    oriC <- 3886105
    terB <- 1661421

    if (c =="oriC") {
        new.origin <- oriC
    } else if (c == 'terB') {
        new.origin <- terB
    } else {
        stop("only oriC or terB are allowed inputs")
    }
    
    if (new.origin >= midpoint) {
        L <- new.origin - midpoint
        ifelse(my.position > L, my.position - new.origin, GENOME.LENGTH - new.origin + my.position)
    } else { ## midpoint is greater than new.origin.
        L <- midpoint + new.origin
        ifelse(my.position > L, my.position - GENOME.LENGTH - new.origin, my.position - new.origin)
    }
}

## look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.
## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length).
calc.cumulative.muts <- function(d, normalization.constant=NA) {

    cumsum.per.pop.helper.func <- function(pop) {
        finalgen <- 6.3 ## this is outside of the data collection
        ## for nice plotting (final generation in mutation.data is 6.275).

        ## This constant is to make sure that all pops are in the levels
        ## of the Population factor after mergers, etc.
        pop.levels <- c("Ara-5","Ara-6", "Ara+1", "Ara+2",
                        "Ara+4", "Ara+5", "Ara-1", "Ara-2",
                        "Ara-3", "Ara-4", "Ara+3", "Ara+6")
        
        df <- d %>% filter(Population==pop)
        if (nrow(df) == 0) { ## if no mutations in this pop.
            almost.done.df <- tibble(Population = factor(pop, levels = pop.levels),
                                     Generation=finalgen,
                                     count=0,
                                     cs=0)
        } else {
            summary.df <- df %>%
                arrange(t0) %>%
                group_by(Population,Generation) %>%
                summarize(count=n()) %>%
                mutate(cs=cumsum(count)) %>%
                ungroup()
            ## if the final generation is not in ret.df,
            ## then add one final row (for nicer plots).
            final.row.df <- tibble(Population=factor(pop, levels = pop.levels),
                                   Generation=finalgen,
                                   count=max(summary.df$count),
                                   cs=max(summary.df$cs))
            
            almost.done.df <- bind_rows(summary.df, final.row.df)
        }
        ## add an row for Generation == 0 (for nicer plots).
        init.row.df <- tibble(
            Population = factor(pop, levels = pop.levels),
            Generation = 0,
            count = 0,
            cs = 0)
        
        ret.df <- bind_rows(init.row.df,almost.done.df)
        return(ret.df)
    }
    
    ## if normalization.constant is not provided, then
    ## calculate based on gene length by default.
    if (is.na(normalization.constant)) {
        my.genes <- d %>% dplyr::select(Gene,gene_length) %>% distinct()
        normalization.constant <- sum(my.genes$gene_length)
    }
    
    c.dat <- map_dfr(.x=levels(d$Population),
                     .f=cumsum.per.pop.helper.func) %>%
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
## Throughout, plots use the minimum subsample size to subsample the null distribution,
## to increase the variance in order to make a conservative comparison.
plot.base.layer <- function(data, subset.size=50, N=1000, alpha = 0.05, normalization.constant=NA, my.color="gray") {

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
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(-alpha/2) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)

    p <- ggplot(filtered.trajectories,aes(x=Generation,y=normalized.cs)) +
        ylab('Cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14))
    return(p)                
}

## add a base layer to a plot. used in Imodulon code.
add.base.layer <- function(p, data, my.color, subset.size=50, N=1000, alpha = 0.05, normalization.constant=NA) {
    
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
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(-alpha/2) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)

    p <- p + geom_point(data=filtered.trajectories,
                        aes(x=Generation, y=normalized.cs),
                        size=0.2, color=my.color)
    return(p)                
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color) {
    p <- p +
        geom_point(data=layer.df,
                   aes(x=Generation,y=normalized.cs),
                   color=my.color, size=0.2) +
        geom_step(data=layer.df, aes(x=Generation,y=normalized.cs),
                  size=0.2, color=my.color)
    return(p)
}

## calculate cumulative numbers of mutations in each category.
## for vanilla plotting, without null distributions, as plotted by
## plot.base.layer.
plot.cumulative.muts <- function(mut.data, my.color="black") {
    p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs)) +
        ylab('Cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)
    return(p)
}

## calculate derivative of cumulative numbers of mutations in each category.
plot.slope.of.cumulative.muts <- function(mut.data, my.color="black") {
    p <- ggplot(mut.data,aes(x=Generation,y=D.normalized.cs)) +
        ylab('Slope of Cumulative number of mutations (normalized)') +
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
add.slope.of.cumulative.mut.layer <- function(p, layer.df, my.color) {
    p <- p +
        geom_point(data=layer.df, aes(x=Generation,y=D.normalized.cs), color=my.color, size=0.2) +
        geom_step(data=layer.df, aes(x=Generation,y=D.normalized.cs), color=my.color, size=0.2) +
        geom_smooth(data=layer.df,size=0.2, color=my.color)
    return(p)
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

## This plot visualizes a two-tailed test (alpha = 0.05)
## against a bootstrapped null distribution for the derivative of cumulative mutations.
## This is what we want to use for publication.
plot.slope.of.base.layer <- function(data, subset.size=300, N=1000, alpha = 0.05, normalization.constant=NA) {

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations,
    ## and then its derivative.
    generate.slope.of.cumulative.mut.subset <- function(idx) {
        rando.genes <- sample(unique(data$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        D.c.mut.subset <- calc.cumulative.muts(mut.subset, normalization.constant) %>%
            calc.slope.of.cumulative.muts() %>%
            mutate(bootstrap_replicate=idx)
        return(D.c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.slope.of.cumulative.mut.subset)

    ## filter out the top alpha/2 and bottom alpha/2 trajectories from each population,
    ## for a two-sided test. default is alpha == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.D.norm.cs=max(D.normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(alpha/2) %>%
        dplyr::select(-final.D.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(-alpha/2) %>%
        dplyr::select(-final.D.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)

    p <- ggplot(filtered.trajectories,aes(x=Generation,y=D.normalized.cs)) +
        ylab('slope of cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color='gray') +
        geom_smooth() +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14))
    return(p)
}

##########################################################################
## DATA ANALYSIS
##########################################################################

## import thetaS estimates.
thetaS.estimates <- read.csv("../data/Martincorena_Maddamsetti_thetaS_estimates.csv")

## get the lengths of all genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
##Do by running:
##python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
    mutate(gene_length=strtoi(gene_length)) %>%
    mutate(oriC_start=rotate.REL606.chr(start,"oriC")) %>%
    mutate(oriC_end=rotate.REL606.chr(end,"oriC")) %>%
    ## join thetaS estimates.
    left_join(thetaS.estimates)

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
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops))) %>%
    ## This is for plotting mutation biases around oriC.
    mutate(oriC.coordinate=rotate.REL606.chr(Position,"oriC")) %>%
    mutate(terB.coordinate=rotate.REL606.chr(Position,"terB"))

## IMPORTANT: DEBUG THIS MERGE-- SHOULD MERGE ON LOCUS_TAG RATHER THAN GENE?
##test <- inner_join(mutation.data,REL606.genes)
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
gene.mutation.data <- gene.mutation.data %>%
    filter(!(Gene %in% duplicate.genes$Gene))

##########################################################################
## parallelism in dS at the same position in the same population
## we get exactly one gene: ydfQ--
## BUT THIS IS A BUG TO BE FIXED, CAUSED BY TWO LOCI WITH THE SAME BLATTNER NUMBER.
bug.to.fix <- gene.dS.mutation.data %>% group_by(Population,Gene,Position) %>% summarize(count=n()) %>% arrange(desc(count)) %>% filter(count>1)
buggy.ydfQ.mutations <- gene.dS.mutation.data %>% filter(Gene=='ydfQ')

##########################################################################################
## MUTATION BIAS ANALYSIS.
## for indels and structural variation, we cannot distinguish between
## mutation hotspots vs. selection.

## EVIDENCE OF STRAND-SPECIFIC BIAS.

## There is evidence of a strand-specific mutation bias on genes
## in the LTEE. Do genes on the lagging strand have a different number of
## mutations compared to the lagging strand?
## Will have to normalize by number of genes in each class.

## The key is that lagging/leading strands flip at the replication origin.
## so the asymmetry on each strand over the origin SHOWS the strand-specfic bias.
## the ratio of total mutations per strand on each side of the origin should give
## an estimate of the strength of this bias.

summed.strand.mut.plot <- ggplot(gene.mutation.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) +
    theme_classic() +
    facet_grid(strand~.,scales="fixed") 
ggsave("../results/figures/summed.strand-mutation-bias-histogram.pdf",summed.strand.mut.plot,width=11,height=8)

araplus3.gene.mutation.data <- gene.mutation.data %>% filter(Population=='Ara+3')
araplus3.strand.mut.plot <- ggplot(araplus3.gene.mutation.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) +
    theme_classic() +
    facet_grid(strand~.,scales="fixed") 
ggsave("../results/figures/ara+3.strand-mutation-bias-histogram.pdf",araplus3.strand.mut.plot,width=11,height=8)

## Variable names assume that the lagging strand has more mutations--
## I can't tell based on arbitrary orientation labeling conventions for
## oriC.coordinate and strand!

leadingstrand.gene.mut.data <- gene.mutation.data %>%
    filter(
    ((oriC.coordinate > 0) & (strand == 1)) | ((oriC.coordinate < 0) & (strand == -1))
    )
laggingstrand.gene.mut.data <- gene.mutation.data %>%
    filter(
    ((oriC.coordinate > 0) & (strand == -1)) | ((oriC.coordinate < 0) & (strand == 1))
    )

nrow(leadingstrand.gene.mut.data)
nrow(laggingstrand.gene.mut.data)

dS.leadingstrand.gene.mut.data <- leadingstrand.gene.mut.data %>%
                                         filter(Annotation=='synonymous')
dS.laggingstrand.gene.mut.data <- laggingstrand.gene.mut.data %>%
                                         filter(Annotation=='synonymous')


araplus3.leadingstrand.gene.mut.data <- leadingstrand.gene.mut.data %>%
    filter(Population=="Ara+3")
araplus3.laggingstrand.gene.mut.data <- laggingstrand.gene.mut.data %>%
    filter(Population=="Ara+3")

dS.araplus3.leadingstrand.gene.mut.data <- araplus3.leadingstrand.gene.mut.data %>%
                                         filter(Annotation=='synonymous')
dS.araplus3.laggingstrand.gene.mut.data <- araplus3.laggingstrand.gene.mut.data %>%
                                         filter(Annotation=='synonymous')

nrow(dS.araplus3.leadingstrand.gene.mut.data)
nrow(dS.araplus3.laggingstrand.gene.mut.data)

leadingstrand.REL606.genes <- REL606.genes %>%
    filter(
    ((oriC_start > 0) & (strand == 1)) | ((oriC_start < 0) & (strand == -1))
    )
laggingstrand.REL606.genes <- REL606.genes %>%
    filter(
    ((oriC_start > 0) & (strand == -1)) | ((oriC_start < 0) & (strand == 1))
    )


leading.gene.target <- sum(leadingstrand.REL606.genes$gene_length)
lagging.gene.target <- sum(laggingstrand.REL606.genes$gene_length)


nrow(leadingstrand.gene.mut.data)/nrow(laggingstrand.gene.mut.data)

nrow(dS.araplus3.leadingstrand.gene.mut.data)/nrow(dS.araplus3.laggingstrand.gene.mut.data)
leading.gene.target/lagging.gene.target

dS.strand.mut.vec <- c(nrow(dS.laggingstrand.gene.mut.data),
                            nrow(dS.leadingstrand.gene.mut.data))

binom.test(x=dS.strand.mut.vec,p=(lagging.gene.target/leading.gene.target))

strand.mut.vec <- c(nrow(laggingstrand.gene.mut.data),
                            nrow(leadingstrand.gene.mut.data))

binom.test(x=strand.mut.vec,p=(lagging.gene.target/leading.gene.target))


## EVIDENCE OF WAVE-PATTERN MUTATION BIAS.

## LOOKS LIKES GENOME-WIDE MUTATION BIAS IN ARA+3, but not others!
## look at mutation density over the chromosome.
## In Ara+3, we see a wave pattern, as reported by Pat Foster's
## group and by Vaughn Cooper's group!
## By looking at mutations in topA, fis, and dusB, Ara+3 has only
## one synonymous mutation in those three genes despite being an
## early mutator. Therefore, perhaps the uniformity
## of dS over LTEE genomes reflect unwinding/loosening of chromosomal
## proteins packing up the DNA.

mut.plot <- ggplot(mutation.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/mutation-bias-histogram.pdf",mut.plot,width=11,height=8)

summed.mut.plot <- ggplot(mutation.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic()
ggsave("../results/figures/summed-mutation-bias-histogram.pdf",mut.plot,width=11,height=8)

## plot point mutations over the genome: dN, dS, nonsense, noncoding
point.mut.data <- mutation.data %>%
     filter(Annotation %in% c('missense', 'synonymous', 'nonsense', 'noncoding'))

point.mut.plot <- ggplot(point.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/point-mut-bias-histogram.pdf",point.mut.plot,width=11,height=8)


## plot non-point mutations over the genome:
nonpoint.mut.data <- mutation.data %>%
     filter(Annotation %in% c('indel', 'sv'))

nonpoint.mut.plot <- ggplot(nonpoint.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/nonpoint-mut-bias-histogram.pdf",nonpoint.mut.plot,width=11,height=8)

summed.nonpoint.mut.plot <- ggplot(nonpoint.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic()
ggsave("../results/figures/summed-nonpoint-mut-bias-histogram.pdf",summed.nonpoint.mut.plot,width=11,height=8)

## plot the distribution of indels and SV over the genome, separately for the
## different LTEE populations, and altogether.

indel.mut.data <- mutation.data %>%
     filter(Annotation %in% c('indel'))

indel.mut.plot <- ggplot(indel.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/indel-mut-bias-histogram.pdf",indel.mut.plot,width=11,height=8)

summed.indel.mut.plot <- ggplot(indel.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic()
ggsave("../results/figures/summed-indel-mut-bias-histogram.pdf",summed.indel.mut.plot,width=11,height=8)

sv.mut.data <- mutation.data %>%
    filter(Annotation %in% c('sv'))

sv.mut.plot <- ggplot(sv.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/sv-mut-bias-histogram.pdf",sv.mut.plot,width=11,height=8)

summed.sv.mut.plot <- ggplot(sv.mut.data,aes(x=oriC.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic()
ggsave("../results/figures/summed-sv-mut-bias-histogram.pdf",summed.sv.mut.plot,width=11,height=8)

## plot gene length against location in oriC coordinates.
gene.length.location.plot <- ggplot(REL606.genes, aes(x=oriC_start,y=gene_length)) +
    geom_point(size=0.5) + geom_smooth() + theme_classic()
gene.length.location.plot

## Figure 1 will be a conceptual figure of the different evolutionary models
## for different sets of genes.
## (purifying selection, coupon collecting, neutral or mutation accumulation, etc.)

## Figure 2 TODO: 
## combine these results to show distribution of different classes of mutations
## over the genome.

#################################################################################
## Control analysis 1:
## look at the accumulation of stars over time for top genes in the
## Tenaillon et al. (2016) genomics data.
## split data into before 50K and after 50K,
## to ask whether we see continued fine-tuning in these genes, overall.

## base plots of null distribution for comparison.
pre50K.rando.plot <- plot.base.layer(filter(gene.mutation.data,Generation <= 5))
post50K.rando.plot <- plot.base.layer(filter(gene.mutation.data,Generation > 5))

## 1) plot top genes in non-mutators.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv')
top.nonmut.genomics <- top_n(nonmut.genomics, 50, wt=G.score)

top.nonmut.mutation.data <- gene.mutation.data %>%
    filter(Gene %in% top.nonmut.genomics$Gene.name)

pre50K.top.nonmut.data <- top.nonmut.mutation.data %>%
    filter(Generation <= 5)

post50K.top.nonmut.data <- top.nonmut.mutation.data %>%
    filter(Generation > 5)

c.pre50K.top.nonmuts <- calc.cumulative.muts(pre50K.top.nonmut.data) %>%
    filter(Generation<=5)
c.post50K.top.nonmuts <- calc.cumulative.muts(post50K.top.nonmut.data) %>%
    filter(Generation > 5)

S1Fig <- pre50K.rando.plot %>%
    add.cumulative.mut.layer(c.pre50K.top.nonmuts,my.color="black")
ggsave("../results/figures/S1Figure.pdf",S1Fig)
S2Fig <- post50K.rando.plot %>%
    add.cumulative.mut.layer(c.post50K.top.nonmuts,my.color="black")
ggsave("../results/figures/S2Figure.pdf",S2Fig)

## data favors coupon-collecting/mutation accumulation:
## genes under selection before 50K don't look so special after 50K.

## 3) plot top genes in hypermutators.
hypermut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv')
top.hypermut.genomics <- top_n(hypermut.genomics, 50, wt=G.score)

top.hypermut.data <- gene.mutation.data %>%
    filter(Gene %in% top.hypermut.genomics$Gene.name)

pre50K.top.hypermut <- top.hypermut.data %>%
    filter(Generation <= 5)

post50K.top.hypermut <- top.hypermut.data %>%
    filter(Generation > 5)

c.pre50K.top.hypermut <- calc.cumulative.muts(pre50K.top.hypermut)
c.post50K.top.hypermut <- calc.cumulative.muts(post50K.top.hypermut)

S3Fig <- pre50K.rando.plot %>%
    add.cumulative.mut.layer(c.pre50K.top.hypermut,my.color="black")
ggsave("../results/figures/S3Figure.pdf",S3Fig)

S4Fig <- post50K.rando.plot %>%
    add.cumulative.mut.layer(c.post50K.top.hypermut,my.color="black")
ggsave("../results/figures/S4Figure.pdf",S4Fig)

#########################################################################
## PURIFYING SELECTION ANALYSIS.

## Find all genes that have no mutations whatsoever in any population.

## IMPORTANT: most of these genes are very short. Many of them are probably
## depleted by chance, and not are significant after FDR-correction.
## nonetheless, some of them are probably under purifying selection.

## Cross-check with the LTEE Genomics data-- synonymous, intergenic (not in coding region)
## and amplifications are allowed, but no other kinds of mutations.

## I downloaded mutations from https://barricklab.org/shiny/LTEE-Ecoli/,
## skipping SNP synonymous, SNP intergenic, and large amplication mutations.
LTEE.genomics.muts <- read.csv("../data/LTEE-Ecoli-data 2020-02-24 14_35_41.csv") %>%
    ## ignore MOB, DEL, etc. that are annotated as intergenic.
    ## this is to prevent false negatives, i.e. removing genes in the gene_list
    ## that weren't directly affected by the intergenic mutation.
    filter(!(str_detect(html_mutation_annotation,"intergenic")))
    
## remove genes from 'no.mutation.genes' that have mutations in LTEE.genomics.muts.
## make a big string of all mutated genes in the LTEE-genomics data, and search this for
## matches.
LTEE.genomics.mutated.genestr <- str_c(unique(LTEE.genomics.muts$gene_list),collapse = ",")

total.gene.mut.count <- nrow(gene.mutation.data)

no.mutation.genes <- REL606.genes %>%
    ## no hits allowed in the metagenomics.
    filter(!(Gene %in% mutation.data$Gene)) %>%
    ## and no hits (other than dS, amps, and intergenic) allowed in the genomics.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length)) %>%
    ## estimate the probability that these genes were not hit by chance.
    mutate(pval = probability.that.not.hit(gene_length,total.gene.mut.count)) %>%
    mutate(fdr.qval = p.adjust(pval,"fdr"))

write.csv(no.mutation.genes,file="../results/no-mutation-genes.csv")

## Look at genes that only have dS.
### MUTATION DENSITIES SUMMED OVER ALL POPULATIONS.
## IMPORTANT: these densities are summed over ALL LTEE populations,
## so they don't correspond to single LTEE populations.

all.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data,
    c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense")) %>%
    rename(all.mut.count = mut.count) %>%
    rename(all.mut.density = density)

all.except.dS.density <- calc.gene.mutation.density(
    gene.mutation.data,c("sv", "indel", "nonsense", "missense")) %>%
    rename(all.except.dS.mut.count = mut.count) %>%
    rename(all.except.dS.mut.density = density)

no.dS.count <- sum(all.except.dS.density$all.except.dS.mut.count)

## combine these into one dataframe.
gene.mutation.densities <- REL606.genes %>%
    full_join(all.mutation.density) %>%
    full_join(all.except.dS.density)

#### CRITICAL STEP: replace NAs with zeros.
#### We need to keep track of genes that haven't been hit by any mutations
#### in a given mutation class (sv, indels, dN, etc.)
gene.mutation.densities[is.na(gene.mutation.densities)] <- 0
gene.mutation.densities <- tbl_df(gene.mutation.densities)

## genes that are only affected by synonymous mutations.
only.dS.allowed.genes <- gene.mutation.densities %>%
    filter(all.except.dS.mut.count == 0) %>%
    ## no hits (other than dS, amps, and intergenic) allowed in the LTEE genomics
    ## to remove false positives.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length)) %>%
    ## estimate the probability that these genes were only hit by dS by chance.
    mutate(pval = probability.that.not.hit(gene_length,no.dS.count)) %>%
    mutate(fdr.qval = p.adjust(pval,"fdr"))

write.csv(only.dS.allowed.genes,file="../results/only-dS-allowed-genes.csv")

## 105 genes just have dS.
just.dS.genes <- only.dS.allowed.genes %>% filter(!(Gene %in% no.mutation.genes$Gene))

## 65 genes don't have even dS in the metagenomic data.
no.dS.genes <- only.dS.allowed.genes %>% filter(Gene %in% no.mutation.genes$Gene)

## If these sets are under purifying selection, then they should be more essential.
## let's examine essentiality from the KEIO collection.
KEIO.data <- read.csv("../data/KEIO_Essentiality.csv", header=TRUE,as.is=TRUE) %>%
    dplyr::select(-JW_id)

KEIO.gene.mutation.densities <- left_join(gene.mutation.densities,KEIO.data)

purifying1 <- KEIO.gene.mutation.densities %>%
    mutate(maybe.purifying = (Gene %in% only.dS.allowed.genes$Gene))

purifying1.plot <- ggplot(purifying1,aes(x=maybe.purifying,y=Score)) +
    theme_classic() + geom_boxplot()
purifying1.plot

purifying2.plot <- ggplot(purifying1,aes(x=Score,fill=maybe.purifying)) + theme_classic() +
    facet_wrap(.~maybe.purifying) + geom_histogram(bins=10) + guides(fill=FALSE)
purifying2.plot

Fig5 <- plot_grid(purifying1.plot,purifying2.plot,nrow=2)
ggsave("../results/figures/Figure5.pdf",Fig5)

## potentially purifying genes have a higher KEIO essentially score, as we would hope.
pur1 <- filter(purifying1,maybe.purifying==TRUE)
notpur1 <- filter(purifying1,maybe.purifying==FALSE)

wilcox.test(x=pur1$Score,notpur1$Score)

## TODO: group together genes-- are they significant when considered as one
## big mutational target? Look at cliques in figure 5, and the proteins
## annotated as hypothetical proteins/toxin-antitoxins as well.

##########################################################################

## quickly look at mutations in phage in REL606 with orthologs
## reported in Louis-Marie Bobay paper,
##Pervasive domestication of defective prophages by bacteria.
## These positions come from S1 Dataset from that paper, for phage
## that are annotated as coming from REL606.

## label genes and mutations by the phage.

phage1 <- gene.mutation.data %>% filter(start>=536929) %>% filter(end<=559782) %>%
    mutate(phage='phage1')
phage2 <- gene.mutation.data %>% filter(start>=787952) %>% filter(end<=799275) %>%
    mutate(phage='phage2')
## no mutations in phage 3: probably masked.
phage3 <- gene.mutation.data %>% filter(start>=880691) %>% filter(end<=904605)
phage4 <- gene.mutation.data %>% filter(start>=1409170) %>% filter(end<=1425863) %>%
    mutate(phage='phage4')
phage5 <- gene.mutation.data %>% filter(start>=1600805) %>% filter(end<=1628589) %>%
    mutate(phage='phage5')
phage6 <- gene.mutation.data %>% filter(start>=2100409) %>% filter(end<=2122347) %>%
    mutate(phage='phage6')

## look at the genes in REL606. do any look like under purifying selection?
phage1.genes <- REL606.genes %>% filter(start>=536929) %>% filter(end<=559782) %>%
    mutate(phage='phage1')
phage2.genes <- REL606.genes %>% filter(start>=787952) %>% filter(end<=799275) %>%
    mutate(phage='phage2')
phage4.genes <- REL606.genes %>% filter(start>=1409170) %>% filter(end<=1425863) %>%
    mutate(phage='phage4')
phage5.genes <- REL606.genes %>% filter(start>=1600805) %>% filter(end<=1628589) %>%
    mutate(phage='phage5')
phage6.genes <- REL606.genes %>% filter(start>=2100409) %>% filter(end<=2122347) %>%
    mutate(phage='phage6')

phage.subset.size <- min(nrow(phage1.genes),
                    nrow(phage2.genes),
                    nrow(phage4.genes),
                    nrow(phage5.genes),
                    nrow(phage6.genes))

phage.rando.layer <- plot.base.layer(gene.mutation.data,subset.size=phage.subset.size)

phage.genes <- rbind(phage1.genes,phage2.genes,phage4.genes,phage5.genes,phage6.genes)
phage.gene.mut.data <- rbind(phage1,phage2,phage4,phage5,phage6)

putative.phage.purifying.selection <- phage.genes %>%
    filter(!(locus_tag %in% phage.gene.mut.data$locus_tag)) %>%
    mutate(in.no.mut.list=(locus_tag %in% no.dS.genes$locus_tag))

no.dS.phage.genes <- phage.genes %>% filter(locus_tag %in% no.dS.genes$locus_tag)
only.dS.allowed.phage.genes <- phage.genes %>% filter(locus_tag %in% only.dS.allowed.genes$locus_tag)

c.phage1 <- calc.cumulative.muts(phage1)
c.phage2 <- calc.cumulative.muts(phage2)
c.phage4 <- calc.cumulative.muts(phage4)
c.phage5 <- calc.cumulative.muts(phage5)
c.phage6 <- calc.cumulative.muts(phage6)

c.all.phage <- calc.cumulative.muts(phage.gene.mut.data)

## plot for phage.
## interesting, but will have to do stats more carefully.
phage.plot <- phage.rando.layer %>%
    ## the black layer is all the other ones together.
    add.cumulative.mut.layer(c.all.phage,my.color="black") %>%
    add.cumulative.mut.layer(c.phage1,my.color="red") %>%
    add.cumulative.mut.layer(c.phage2,my.color="blue") %>%
    add.cumulative.mut.layer(c.phage4,my.color="green") %>%
    add.cumulative.mut.layer(c.phage5,my.color="yellow") %>%
    add.cumulative.mut.layer(c.phage6,my.color="orange")
ggsave(phage.plot,filename='../results/figures/phage-plot.pdf')

##########################################################################
## look at accumulation of stars over time for genes in different transcriptional
## modules inferred by Sastry et al. (2020) paper from Bernhard Palsson's group.

## Questions:
## Q1) are the regulators of the I-modulons under stronger selection than
## the genes that they regulate?

## Two ways of testing:
## A) group all regulators together, against all regulated genes.
## B) paired comparison. compare regulator against regulated genes for each I-modulon.

## Q2) Do any I-modulons show evidence of positive or purifying selection?
## I will have to do some kind of FDR correction for multiple hypothesis testing.

## I made this file by hand, by going through the I-modulons in imodulon_gene_names.txt,
## and the regulator annotations in modulon.pdf, both in precise-db-repo.
Imodulons.to.regulators <- read.csv("../data/rohans-I-modulons-to-regulators.csv")

Imodulon.regulators <-Imodulons.to.regulators %>%
    filter(!(is.na(regulator)))
    
Imodulon.regulator.mut.data <- gene.mutation.data %>%
    filter(Gene %in% Imodulon.regulators$regulator)

c.Imodulon.regulators <- calc.cumulative.muts(Imodulon.regulator.mut.data)

## calculate more rigorous statistics than the figures.
## IMPORTANT BUG: WHY ARE ONLY 7 OF 12 POPS IN the PVAL DATAFRAME ???
Imodulon.regulator.pvals <- calculate.trajectory.tail.probs(gene.mutation.data, unique(Imodulon.regulators$regulator))

## Now look at genes that are regulated within Imodulons.
## I expect relaxed or purifying selection overall.

## this file was generated by reformat-I-modulons.py
genes.to.Imodulons <- read.csv("../results/genes-to-I-modulons.csv")

Imodulon.regulated <- genes.to.Imodulons %>% filter(!(is.na(Gene))) %>%
    ## remove any I-modulon regulators from the I-modulon genes to
    ## make an orthogonal comparison.
    filter(!(Gene %in% Imodulon.regulators$regulator))
    
Imodulon.regulated.mut.data <- gene.mutation.data %>%
    filter(Gene %in% Imodulon.regulated$Gene)

c.Imodulon.regulated <- calc.cumulative.muts(Imodulon.regulated.mut.data)

Imodulon.regulators.base.layer <- plot.base.layer(
    gene.mutation.data,
    subset.size=length(unique(Imodulon.regulators$regulator)))

## Figure for paper:  compare I-modulon regulators to the genes they regulate.
Imodulon.plot <- Imodulon.regulators.base.layer %>% ## null for regulators
    add.base.layer(gene.mutation.data, ## add null for regulated genes
                   subset.size=length(unique(Imodulon.regulated$Gene)),
                   my.color="pink") %>%
    add.cumulative.mut.layer(c.Imodulon.regulators, my.color="black") %>%
    add.cumulative.mut.layer(c.Imodulon.regulated, my.color="red")
ggsave("../results/figures/Imodulon-plot.pdf",Imodulon.plot)


## Make plots for each I-modulon.

## This helper function refers to several global variables.
## TODO: Best to wrap this into the context of a larger function later,
## to maintain modularity.
make.modulon.plots.helper <- function(my.I.modulon) {
    
    my.modulon.name <- unique(my.I.modulon$I.modulon)
    modulon.text <- paste(my.modulon.name,"I-modulon")
    print(modulon.text) ## to help with debugging.

    my.regulators <- my.I.modulon %>% filter(!(is.na(regulator)))
    regulator.size <- length(unique(my.regulators$regulator))
    
    my.modulon.genes <- genes.to.Imodulons %>%
        filter(I.modulon == my.modulon.name) %>%
        filter(!(is.na(Gene))) %>%
        ## remove any regulators of this I-modulon (if they exist)
        ## to make an orthogonal comparison.
        filter(!(Gene %in% my.regulators$regulator))
    
    my.modulon.mut.data <- gene.mutation.data %>%
        filter(Gene %in% my.modulon.genes$Gene)

    c.my.modulon.muts <- calc.cumulative.muts(my.modulon.mut.data)
    ## for the plots, subsample based on the cardinality of the I-modulon.
    modulon.size <- length(unique(my.modulon.genes$Gene))
    modulon.base.layer <- plot.base.layer(gene.mutation.data,
                                          subset.size=modulon.size,
                                          my.color="pink")
    p <- modulon.base.layer %>%
        add.cumulative.mut.layer(c.my.modulon.muts,my.color="red")

    ## if any regulators exist, add layers for mutations in those genes.
    if (regulator.size > 0) {
        my.regulator.mut.data <- gene.mutation.data %>%
            filter(Gene %in% my.regulators$regulator)
        c.my.regulators <- calc.cumulative.muts(my.regulator.mut.data)

        p <- p %>% 
            add.base.layer(gene.mutation.data,my.color="grey",
                           subset.size=regulator.size) %>%
            add.cumulative.mut.layer(c.my.regulators,my.color='black')
    }
    ## add a title to the plots.
    p <- p + ggtitle(modulon.text)
    print(p)
}

## too big for PDF! plot to separate jpegs.
jpeg("../results/figures/I-modulon-plots/plot-%d.jpeg")
## split into groups by I-modulon, and map them to the plotting helper.
Imodulons.to.regulators %>% group_split(I.modulon) %>%
    map(.f=make.modulon.plots.helper)
dev.off()

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

proteome.subset.size <- min(length(unique(A.sector.mut.data$Gene)),
                            length(unique(S.sector.mut.data$Gene)),
                            length(unique(O.sector.mut.data$Gene)),
                            length(unique(U.sector.mut.data$Gene)),
                            length(unique(R.sector.mut.data$Gene)),
                            length(unique(C.sector.mut.data$Gene)))

proteome.rando.layer <- plot.base.layer(gene.mutation.data,subset.size=proteome.subset.size)


c.A.muts <- calc.cumulative.muts(A.sector.mut.data)
c.S.muts <- calc.cumulative.muts(S.sector.mut.data)
c.O.muts <- calc.cumulative.muts(O.sector.mut.data)
c.U.muts <- calc.cumulative.muts(U.sector.mut.data)
c.R.muts <- calc.cumulative.muts(R.sector.mut.data)
c.C.muts <- calc.cumulative.muts(C.sector.mut.data)

## plot for proteome sectors.
sector.plot <- proteome.rando.layer %>%
    add.cumulative.mut.layer(c.A.muts, my.color="black") %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red") %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue") %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green") %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow") %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange")
ggsave(sector.plot,filename='../results/figures/sector-plot.pdf')

##########################################################################
## look at accumulation of stars over time for genes in different eigengenes
## inferred by Wytock and Motter (2018).

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

eigen.subset.size <- min(length(unique(eigengene1.mut.data$Gene)),
                         length(unique(eigengene2.mut.data$Gene)),
                         length(unique(eigengene3.mut.data$Gene)),
                         length(unique(eigengene4.mut.data$Gene)),
                         length(unique(eigengene5.mut.data$Gene)),
                         length(unique(eigengene6.mut.data$Gene)),
                         length(unique(eigengene7.mut.data$Gene)),
                         length(unique(eigengene8.mut.data$Gene)),
                         length(unique(eigengene9.mut.data$Gene)))

eigen.rando.layer <- plot.base.layer(gene.mutation.data,subset.size=eigen.subset.size)

c.eigen1.muts <- calc.cumulative.muts(eigengene1.mut.data)
c.eigen2.muts <- calc.cumulative.muts(eigengene2.mut.data)
c.eigen3.muts <- calc.cumulative.muts(eigengene3.mut.data)
c.eigen4.muts <- calc.cumulative.muts(eigengene4.mut.data)
c.eigen5.muts <- calc.cumulative.muts(eigengene5.mut.data)
c.eigen6.muts <- calc.cumulative.muts(eigengene6.mut.data)
c.eigen7.muts <- calc.cumulative.muts(eigengene7.mut.data)
c.eigen8.muts <- calc.cumulative.muts(eigengene8.mut.data)
c.eigen9.muts <- calc.cumulative.muts(eigengene9.mut.data)

eigen.plot <- eigen.rando.layer %>%
    add.cumulative.mut.layer(c.eigen1.muts, my.color="red") %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange") %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow") %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green") %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan") %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue") %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet") %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink") %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black")
ggsave(eigen.plot,filename='../results/figures/eigen-plot.pdf')

## why does Ara+3 an upswing of mutations in eigengene 6?
## Seems to be entirely driven by an excess of mutations in entF.
## Looks like a bona fide example of historical contingency!
eigengene6.mut.data %>% filter(Population=='Ara+3') %>%
    group_by(Gene) %>% summarize(count=n())

##########################################################################################
## SUPPLEMENTARY INFORMATION
##########################################################################################

## reanalysis of synonymous variation in natural populations.
## looks like problems with the K-S test set up?

## investigate dS (and other classes of mutations) across the genome
## in the metagenomics data.
## revamp code from my 2015 Mol. Biol. Evol. paper.
## in short, cannot reject null model that dS is uniform over the genome.
## and dN fits the null extremely well-- even better than dS!
## Probably because there are 3 times as many dN as dS throughout the
## experiment... but could the mutation bias seen in Ara+3 be a factor?

## IMPORTANT ISSUE: changing the CDF (by changing rank on x-axis),
## while keeping the same set of probability masses for each gene
## changes K-S test statistics! How do I interpret this?
## Perhaps because K-S is for continuous and not discrete distributions?

########################################################################
## DO NOT replace this with thetaS analysis-- that will only work for core genes!
ks.analysis <- function(the.data, REL606.genes, order_by_oriC=FALSE) {
    ## For each set of data (all data, non-mutators, MMR mutators, mutT mutators)
    ## do the following: 1) make a uniform cdf on mutation rate per base,
    ## 2) make an empirical cdf of mutations per gene.
    ## do K-S tests for goodness of fit of the empirical cdf with
    ## the uniform cdf.
    
    hit.genes.df <- the.data %>%
        group_by(locus_tag, Gene, gene_length, oriC_start) %>%
        summarize(hits=n()) %>%
        ungroup()
    
    ## have to do it this way, so that zeros are included.
    hit.genes.df <- full_join(REL606.genes,hit.genes.df) %>% replace_na(list(hits=0)) %>%
        arrange(desc(gene_length))

    if (order_by_oriC) { ## oriC will be roughly in the middle
        hit.genes.df <- hit.genes.df %>% arrange(oriC_start)
    }

    
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
                                  null=null.cdf,
                                  oriC_start=hit.genes.df$oriC_start)
    return(results.to.plot)
}

make.KS.Figure <- function(the.results.to.plot, order_by_oriC=FALSE) {
    
    if (order_by_oriC) {
        the.results.to.plot <- the.results.to.plot %>% arrange(oriC_start)
    } else {
        the.results.to.plot <- the.results.to.plot %>% arrange(gene_length)
    }
    
    ## for plotting convenience, add an index to the data frame.
    the.results.to.plot$index <- seq_len(nrow(the.results.to.plot))
  
    p <- ggplot(the.results.to.plot, aes(x=index)) +
        geom_line(aes(y=empirical), colour="red") + 
        geom_line(aes(y=null), linetype=2) + 
        scale_y_continuous('Cumulative proportion of mutations',limits=c(0,1)) +
        theme_classic() +
        theme(axis.title=element_text(size=18),axis.text=element_text(size=12))

    if (order_by_oriC) {
        p <- p + scale_x_continuous('Genes ranked by chromosomal location',
                                    limits=c(0,4400))
    } else {
        p <- p + scale_x_continuous('Genes ranked by length',limits=c(0,4400))
    }
    return(p)
}

## IMPORTANT NOTE: THIS WILL ONLY WORK ON CORE GENES.
## TODO: add use.maddamsetti parameter to use one or other set of parameter estimates.
## for experimenting with how the ranking on the x-axis changes the statistics,
## I added the rev parameter to reverse the order of genes on the x-axis.
thetaS.KS.analysis <- function(the.data, REL606.genes, rank_by="length") {

    ## rank_by can have three values: "length", "thetaS", or "oriC".
    stopifnot(rank_by %in% c("length","thetaS","oriC"))
    
    ## 1) make an empirical cdf of mutations per core gene.
    ## do K-S tests for goodness of fit of the empirical cdf with cdfs for
    ## thetaS.
    
    hit.genes.df <- the.data %>%
        group_by(locus_tag, Gene, gene_length, oriC_start) %>%
        summarize(hits=n()) %>%
        ungroup()

    ## have to do it this way, so that zeros are included.
    hit.genes.df <- full_join(REL606.genes,hit.genes.df) %>%
        mutate(thetaS=Martincorena_thetaS) %>%
        ##mutate(thetaS=Maddamsetti_thetaS) %>%
        filter(!(is.na(thetaS))) %>% ## only keep core genes.
        replace_na(list(hits=0))

    if (rank_by == "oriC") {
        hit.genes.df <- hit.genes.df %>%
            arrange(oriC_start) ## oriC will be roughly in the middle
    } else if (rank_by == "length") {
        hit.genes.df <- hit.genes.df %>%
            arrange(gene_length)
    } else if (rank_by == "thetaS") {
        hit.genes.df <- hit.genes.df %>%
            arrange(thetaS)
    } else {
        stop("error in thetaS.KS.analysis.")
    }
    
    ## Calculate the empirical distribution of synonymous substitutions per gene.
    hit.genes.length <- sum(hit.genes.df$gene_length)
    mutation.total <- sum(hit.genes.df$hits)
    empirical.cdf <- cumsum(hit.genes.df$hits)/mutation.total
    ## Null hypothesis: probability of a mutation per base is uniform.
    null.cdf <- cumsum(hit.genes.df$gene_length)/hit.genes.length
    ## Alternative hypothesis: mutation rate is proportional to thetaS.
    ## alternative 1: thetaS is the mutation rate per base pair.
    thetaS.is.per.bp.total <- sum(hit.genes.df$thetaS * hit.genes.df$gene_length)
    alt1.cdf <- cumsum(hit.genes.df$thetaS * hit.genes.df$gene_length /thetaS.is.per.bp.total)
    ## alternative 2: thetaS is the mutation rate per gene.
    thetaS.is.per.gene.total <- sum(hit.genes.df$thetaS)
    alt2.cdf <- cumsum(hit.genes.df$thetaS/thetaS.is.per.gene.total)
    
    ## Do Kolmogorov-Smirnov tests for goodness of fit.
    print("uniform mutation rate hypothesis")
    print(ks.test(empirical.cdf, null.cdf, simulate.p.value=TRUE))
    print("mutation rate per bp is proportional to thetaS hypothesis")
    print(ks.test(empirical.cdf, alt1.cdf, simulate.p.value=TRUE))
    print("mutation rate per gene is proportional to thetaS hypothesis")
    print(ks.test(empirical.cdf, alt2.cdf, simulate.p.value=TRUE))
    
    results.to.plot <- data.frame(locus_tag=hit.genes.df$locus_tag,
                                  Gene=hit.genes.df$Gene,
                                  gene_length=hit.genes.df$gene_length,
                                  thetaS=hit.genes.df$thetaS,
                                  empirical=empirical.cdf,
                                  null=null.cdf,
                                  alt1=alt1.cdf,
                                  alt2=alt2.cdf,
                                  oriC_start=hit.genes.df$oriC_start)
    return(results.to.plot)
}

## for experimenting with how the ranking on the x-axis changes the statistics,
## I added the rev parameter to reverse the order of genes on the x-axis.
make.thetaS.KS.Figure <- function(the.results.to.plot, rank_by="length") {

    ## rank_by can have three values: "length", "thetaS", or "oriC".
    stopifnot(rank_by %in% c("length","thetaS","oriC"))

    if (rank_by == "oriC") {
        ## oriC will be roughly in the middle
        the.results.to.plot <- the.results.to.plot %>% arrange(oriC_start)
    } else if (rank_by == "length") {
        the.results.to.plot <- the.results.to.plot %>% arrange(gene_length)
    } else if (rank_by == "thetaS") {
        the.results.to.plot <- the.results.to.plot %>% arrange(thetaS)
    } else {
        stop("error in make.KS.Figure.")
    }

    ## for plotting convenience, add an index to the data frame.
    the.results.to.plot$index <- seq_len(nrow(the.results.to.plot))
    
    plot <- ggplot(the.results.to.plot, aes(x=index)) +
        geom_line(aes(y=empirical), colour="red") + 
        geom_line(aes(y=null), linetype=2) +
        geom_line(aes(y=alt1), linetype='dotted') +
        geom_line(aes(y=alt2), linetype='dotted') + 
        scale_y_continuous('Cumulative proportion of mutations',limits=c(0,1)) +
        theme_classic() +
        theme(axis.title=element_text(size=18),axis.text=element_text(size=12))

    if (rank_by == "oriC") {
        ## find the index for gidA and mioC, which sandwich oriC.
        ## use these to plot the location of oriC.
        gidA.index <- filter(the.results.to.plot,Gene=='gidA')$index
        mioC.index <- filter(the.results.to.plot,Gene=='mioC')$index
        plot <- plot + scale_x_continuous('Genes ranked by chromosomal location') + geom_vline(xintercept=(gidA.index+mioC.index)/2,linetype='dotted',color='grey')
    } else if (rank_by == "length") {
        plot <- plot + scale_x_continuous('Genes ranked by length')
    } else if (rank_by == "thetaS") {
        plot <- plot + scale_x_continuous('Genes ranked by thetaS')
    } else {
        stop("error 2 in make.KS.Figure.")
    }
    
    return(plot)
}




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

## we treat nonsense as a nonsynonymous mutation-- for comparison to Tenaillon 2016 data.
gene.nodNdS.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation %in% c("sv", "indel", "noncoding"))





## examine all mutations over genes in the genome.
cumsum.all.over.metagenome <- ks.analysis(gene.only.mutation.data,REL606.genes)

## NOTE: thetaS KS analysis will only work for core genes!
## IMPORTANT!! RESULTS DEPEND ON X-AXIS ORDERING!
cumsum.dS.core1 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"length")
cumsum.dS.core2 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"oriC")
cumsum.dS.core3 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"thetaS")

cumsum.dN.core1 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"length")
cumsum.dN.core2 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"oriC")
cumsum.dN.core3 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"thetaS")

dS.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dS.core1,"length")
ggsave("../results/figures/dS_thetaS_plot1.pdf",dS.thetaS.plot1)
dS.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dS.core2,"oriC")
ggsave("../results/figures/dS_thetaS_plot2.pdf",dS.thetaS.plot2)
dS.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dS.core3,"thetaS")
ggsave("../results/figures/dS_thetaS_plot3.pdf",dS.thetaS.plot3)

dN.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dN.core1,"length")
ggsave("../results/figures/dN_thetaS_plot1.pdf",dN.thetaS.plot1)
dN.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dN.core2,"oriC")
ggsave("../results/figures/dN_thetaS_plot2.pdf",dN.thetaS.plot2)
dN.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dN.core3,"thetaS")
ggsave("../results/figures/dN_thetaS_plot3.pdf",dN.thetaS.plot3)

##########################

## not significant: but marginal result: p-value = 0.095.
cumsum.all.over.metagenome.by.oriC <- ks.analysis(gene.only.mutation.data,REL606.genes,TRUE)
all.KS.plot.by.oriC <- make.KS.Figure(cumsum.all.over.metagenome.by.oriC, TRUE)
ggsave("../results/figures/all_KS.plot_by_oriC.pdf",all.KS.plot.by.oriC)

## plot all mutations for Ara+3.
ara.plus3.all.muts <- filter(gene.only.mutation.data,Population=='Ara+3')
cumsum.ara.plus3.core.metagenome <- thetaS.KS.analysis(ara.plus3.all.muts,REL606.genes)
ara.plus3.all.thetaS.KS.plot <- make.thetaS.KS.Figure(cumsum.ara.plus3.core.metagenome)

cumsum.all.over.ara.plus3 <- ks.analysis(ara.plus3.all.muts,REL606.genes)
ara.plus3.all.KS.plot <- make.KS.Figure(cumsum.all.over.ara.plus3)
##ggsave("../results/figures/Ara+3_KS.plot.pdf",ara.plus3.all.KS.plot)

cumsum.all.over.ara.plus3.by.oriC <- ks.analysis(ara.plus3.all.muts,REL606.genes,TRUE)
ara.plus3.all.KS.plot.by.oriC <- make.KS.Figure(cumsum.all.over.ara.plus3.by.oriC, TRUE)
##ggsave("../results/figures/Ara+3_KS.plot.pdf",ara.plus3.all.KS.plot)

## plot all mutations for Ara+6.
ara.plus6.all.muts <- filter(gene.only.mutation.data,Population=='Ara+6')

cumsum.all.over.ara.plus6 <- ks.analysis(ara.plus6.all.muts,REL606.genes)
ara.plus6.all.KS.plot <- make.KS.Figure(cumsum.all.over.ara.plus6)
##ggsave("../results/figures/Ara+3_KS.plot.pdf",ara.plus3.all.KS.plot)

cumsum.all.over.ara.plus6.by.oriC <- ks.analysis(ara.plus6.all.muts,REL606.genes,TRUE)
ara.plus6.all.KS.plot.by.oriC <- make.KS.Figure(cumsum.all.over.ara.plus6.by.oriC, TRUE)
##ggsave("../results/figures/Ara+3_KS.plot.pdf",ara.plus3.all.KS.plot)


## examine dS over the genome.
cumsum.dS.over.metagenome <- ks.analysis(gene.dS.mutation.data,REL606.genes)
dS.KS.plot <- make.KS.Figure(cumsum.dS.over.metagenome)
ggsave("../results/figures/dS_KS.plot.pdf",dS.KS.plot)

## dS is significantly different from null, when ordered by chromosomal location!
## p = 0.002231.
cumsum.dS.over.metagenome.by.oriC <- ks.analysis(gene.dS.mutation.data,REL606.genes,TRUE)
dS.KS.plot.by.oriC <- make.KS.Figure(cumsum.dS.over.metagenome.by.oriC,TRUE)
ggsave("../results/figures/dS_KS.plot_by_oriC.pdf",dS.KS.plot.by.oriC)

## plot dS for Ara+3.
ara.plus3.dS.muts <- filter(gene.dS.mutation.data,Population=='Ara+3')

cumsum.dS.over.ara.plus3 <- ks.analysis(ara.plus3.dS.muts,REL606.genes)
ara.plus3.dS.KS.plot <- make.KS.Figure(cumsum.dS.over.ara.plus3)
ggsave("../results/figures/Ara+3_dS_KS.plot.pdf",ara.plus3.dS.KS.plot)

cumsum.dS.over.ara.plus3.by.oriC <- ks.analysis(ara.plus3.dS.muts,REL606.genes,TRUE)
ara.plus3.dS.KS.plot.by.oriC <- make.KS.Figure(cumsum.dS.over.ara.plus3.by.oriC,TRUE)
ggsave("../results/figures/Ara+3_dS_KS_by_oriC.plot.pdf",ara.plus3.dS.KS.plot.by.oriC)


cumsum.ara.plus3.dS.metagenome <- thetaS.KS.analysis(ara.plus3.dS.muts,REL606.genes)
ara.plus3.dS.thetaS.KS.plot <- make.thetaS.KS.Figure(cumsum.ara.plus3.dS.metagenome)


## examine dN over the genome.
## dN fits the null even better than dS! But note that
## there are 18493 dN in the data, and 6792 dS in the data.
## so the better fit is probably best explained by the larger sample size.
cumsum.dN.over.metagenome <- ks.analysis(gene.dN.mutation.data,REL606.genes)
dN.KS.plot <- make.KS.Figure(cumsum.dN.over.metagenome)
cumsum.dN.over.metagenome.by.oriC <- ks.analysis(gene.dN.mutation.data,REL606.genes,TRUE)
dN.KS.plot.by.oriC <- make.KS.Figure(cumsum.dN.over.metagenome.by.oriC,TRUE)
ggsave("../results/figures/dN_KS_plot.pdf",dN.KS.plot)
ggsave("../results/figures/dN_KS_plot_by_oriC.pdf",dN.KS.plot.by.oriC)

## dN is not different from null, when ordered by chromosomal location.
## p = 0.3705.
cumsum.dN.over.metagenome.by.oriC <- ks.analysis(gene.dN.mutation.data,REL606.genes,TRUE)
dN.KS.plot.by.oriC <- make.KS.Figure(cumsum.dN.over.metagenome.by.oriC, TRUE)

## plots for Ara+3.
ara.plus3.dN.muts <- filter(gene.dN.mutation.data,Population=='Ara+3')
cumsum.dN.over.ara.plus3 <- ks.analysis(ara.plus3.dN.muts,REL606.genes)
ara.plus3.dN.KS.plot <- make.KS.Figure(cumsum.dN.over.ara.plus3)

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
## BUT ACTUALLY OCCUR IN PROMOTER REGIONS! I can do this using the gene_start and
## gene_end information.

## examine structural mutations (IS elements) affecting genes.
## RESULT: longer genes are depleted in IS insertions!! (double check if true.)
cumsum.sv.over.metagenome <- ks.analysis(gene.sv.mutation.data,REL606.genes)
make.KS.Figure(cumsum.sv.over.metagenome)

## examine indels.
## RESULT: longer genes are depleted in indels! (double check if true.)
cumsum.indel.over.metagenome <- ks.analysis(gene.indel.mutation.data,REL606.genes)
make.KS.Figure(cumsum.indel.over.metagenome)

## Therefore, structural mutations and indels are probably what are driving the
## differences that I saw in aerobic and anerobic mutations before.

cumsum.nonsense.sv.indels.over.metagenome <- ks.analysis(gene.nonsense.sv.indels.mutation.data, REL606.genes)
make.KS.Figure(cumsum.nonsense.sv.indels.over.metagenome)
