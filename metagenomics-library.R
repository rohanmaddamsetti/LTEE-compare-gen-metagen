## metagenomics-library.R by Rohan Maddamsetti
## This file is a library of common functions used in both paper 2A and paper 2B.

library(tidyverse)
library(cowplot)

##########################################################################
## FUNCTIONS FOR DATA ANALYSIS
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
