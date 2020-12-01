## metagenomics-library.R by Rohan Maddamsetti
## This file is a library of common functions used in both paper 2A and paper 2B.

library(tidyverse)
library(cowplot)

##########################################################################
## FUNCTIONS FOR DATA ANALYSIS
##########################################################################

## function for plotting better y-axis labels.
## see solution here for nice scientific notation on axes.
## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
fancy_scientific <- function(x) {
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

## calculate the probability that a locus (or set of loci) is not hit by mutations,
## assuming uniform mutation rate.
## l is locus length, and n is the number of mutations in the dataset.
uniform.probability.that.not.hit <- function(l,n) {
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

## look at accumulation of stars over time
## in other words, look at the rates at which the mutations occur over time.
## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length).
## If plot.to.end is TRUE, then add one final row.
calc.cumulative.muts <- function(d, d.metadata, plot.to.end=TRUE) {

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
                summarize(count=n(), .groups = "drop_last") %>%
                mutate(cs=cumsum(count)) %>%
                ungroup()
            ## if the final generation is not in ret.df,
            ## then add one final row (for nicer plots).
            final.row.df <- tibble(Population=factor(pop, levels = pop.levels),
                                   Generation=finalgen,
                                   count=max(summary.df$count),
                                   cs=max(summary.df$cs))
            if (plot.to.end) {
                almost.done.df <- bind_rows(summary.df, final.row.df)
            } else {
                almost.done.df <- summary.df
            }
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

    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my.genes <- d.metadata %>% dplyr::select(Gene,gene_length) %>% distinct()
    normalization.constant <- sum(my.genes$gene_length)
    
    c.dat <- map_dfr(.x=levels(d$Population),
                     .f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    
    return(c.dat)
}

###############################
## the next two functions are for estimating local mutation rates by
## splitting the genome into z bins.

## find.bin is useful for sampling random genes while preserving bin identity,
## that is, only sampling genes near the genes of interest.

bin.mutations <- function(mut.data, z) {
    ## z is the number of bins we are using.
    ## count the number of mutations in each bin, M(z_i).
    ## filter araplus3.mut.data using each adjacent pair of fenceposts,
    ## and count the number of rows to get mutations per bin.
    mutations.by.bin.vec <- rep(0,z)

    GENOME.LENGTH <- 4629812
    c <- GENOME.LENGTH/z ## length of each bin
    
    ## define fenceposts for each bin.
    ## this has z+1 entries.
    fenceposts <- seq(0,z) * c
    
    for (i in 1:z) {
        left <- fenceposts[i]
        right <- fenceposts[i+1]
        bin.data <- araplus3.mut.data %>%
            filter(Position >= left & Position < right)
        bin.mut.count <- nrow(bin.data)
        mutations.by.bin.vec[i] <- bin.mut.count
    }
    
    ## assert that all mutations have been assigned to a bin.
    stopifnot(sum(mutations.by.bin.vec) == nrow(araplus3.mut.data))
    
    return(mutations.by.bin.vec)
}

find.bin <- function(locus.row, z) {
    ## take a 1-row dataframe correponding to a REL606 gene,
    ## and return the bin that it belongs to, given z bins.

    GENOME.LENGTH <- 4629812
    c <- GENOME.LENGTH/z ## length of each bin
    
    ## define right-hand fencepost for each bin. 
    rightfencepost <- seq(1,z) * c
    
    for (i in 1:z) {
        if ((locus.row$start < rightfencepost[i]) &
            (locus.row$end < rightfencepost[i]))
            return(i)
    }
    stopifnot(TRUE) ## we should always return a value in the for loop.
    ## There is an unhandled corner case where the gene could straddle a boundary.
    ## So we need to make sure that this doesn't happen in practice. 
    return(-1)
}

###############################

## calculate the tail probabilities of the true cumulative mutation trajectory
## of a given vector of genes (a 'module'), based on resampling
## random sets of genes. Returns the upper tail of null distribution,
## or P(random trajectory >= the actual trajectory).
## Output: a dataframe with three columns: Population, count, p.val
calc.traj.pvals <- function(data, REL606.genes, gene.vec, N=10000, sample.genes.by.location=FALSE) {

    ## check type for gene.vec (e.g., if factor, change to vanilla vector of strings)
    gene.vec <- as.character(gene.vec)
    
    ## each sample has the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)

    if (sample.genes.by.location) {
        ## then sample genes near the genes in the module of interest,
        ## i.e. gene.vec.

        ## for efficiency, pre-calculate gene bin assignments.
        ## use 46 bins so that each bin is ~10000 bp. 
        find.bin.46 <- partial(find.bin,z=46)
        gene.info <- REL606.genes %>%
            select(Gene, locus_tag, blattner, gene_length,
                   product, start, end, strand) %>%
            distinct() ## remove duplicate rows,
        ## as mutations in the same gene have the same gene info
        
        ## output a bin from each row in gene.info.
        bin.list <- gene.info %>%
            split(.$Gene) %>%
            lapply(find.bin.46)
        ## map each gene in the module of interest to their bin.
        bin.vec <- sapply(gene.vec,function(gene) bin.list[[gene]])

        ## associate each gene in the genome to its bin, and add as a column.
        bin.df <- data.frame(Gene=names(bin.list),bin=as.numeric(bin.list))
        
        ## map bins to the genes in those bins (excepting those in gene.vec).
        nested.gene.info.with.bins <- left_join(gene.info, bin.df) %>%
            ## when sampling, exclude the genes in the module of interest.
            filter(!(Gene %in% gene.vec)) %>%
            group_by(bin) %>%
            nest()

        ## add the genes to be sampled as a list column corresponding
        ## to the genes in gene.vec,
        ## such that each gene in the module maps to the other genes in its bin.
        ## then, one of those genes will be sampled to make the random module,
        ## while preserving bin location.
        module.info.with.bins <- left_join(gene.info, bin.df) %>%
            filter(Gene %in% gene.vec) %>%
            left_join(nested.gene.info.with.bins)

        sample.one <- partial(sample_n,size=1)
        
        sample.genes.by.genomebin <- function() {
            ## map each gene in gene.vec to a random gene in its bin
            ## in the genome.
            module.to.random.module <- module.info.with.bins %>%
                mutate(sampled.info = map(data, sample.one)) %>%
                transmute(Gene, sampled.gene = map_chr(sampled.info, function(x) x$Gene))
            return(module.to.random.module$sampled.gene)
        }

        generate.cumulative.mut.subset.by.loc <- function(idx) {
            rando.genes <- sample.genes.by.genomebin()
            mut.subset <- filter(data,Gene %in% rando.genes)
            mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
            c.mut.subset <- calc.cumulative.muts(mut.subset,
                                                 mut.subset.metadata) %>%
                mutate(bootstrap_replicate=idx)
            return(c.mut.subset)
        }

        ## make a dataframe of bootstrapped trajectories.
        ## look at accumulation of stars over time for random subsets of genes.
        bootstrapped.trajectories <- map_dfr(.x=seq_len(N),
                                             .f=generate.cumulative.mut.subset.by.loc)
    } else {
        ## This function takes the index for the current draw, and samples the data,
        ## generating a random gene set for which to calculate cumulative mutations.
        ## IMPORTANT: this function depends on variables defined in
        ## calculate.trajectory.tail.probs.
        generate.cumulative.mut.subset <- function(idx) {
            rando.genes <- base::sample(unique(data$Gene),subset.size)
            mut.subset <- filter(data,Gene %in% rando.genes)
            mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
            c.mut.subset <- calc.cumulative.muts(mut.subset, mut.subset.metadata) %>%
                mutate(bootstrap_replicate=idx)
            return(c.mut.subset)
        }
        
        ## make a dataframe of bootstrapped trajectories.
        ## look at accumulation of stars over time for random subsets of genes.
        bootstrapped.trajectories <- map_dfr(.x=seq_len(N),
                                             .f=generate.cumulative.mut.subset)
    }
    
    gene.vec.data <- data %>% filter(Gene %in% gene.vec)
    data.trajectory <- calc.cumulative.muts(gene.vec.data, gene.vec.metadata)
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
        group_by(Population,.drop=FALSE) %>%
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

## This plot visualizes a two-tailed test (alphaval = 0.05)
## against a bootstrapped null distribution.
## Throughout, plots use the minimum subsample size to subsample the null distribution,
## to increase the variance in order to make a conservative comparison.
plot.base.layer <- function(data, REL606.genes, subset.size=50, N=1000, alphaval = 0.05, my.color="gray") {
    
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## IMPORTANT: this function depends on variables defined in plot.base.layer.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(REL606.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset, mut.subset.metadata) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories from each population,
    ## for a two-sided test. default is alphaval == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
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
              axis.text.y  = element_text(size=14)) +
        scale_y_continuous(labels=fancy_scientific,
                           breaks = scales::extended_breaks(n = 6),
                           limits = c(0, NA))
    return(p)
}

## add a base layer to a plot. used in Imodulon code.
add.base.layer <- function(p, data, REL606.genes, my.color, subset.size=50, N=1000, alphaval = 0.05, normalization.constant=NA) {
    
        ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(REL606.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(REL606.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset, mut.subset.metadata) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.
    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories from each population,
    ## for a two-sided test. default is alphaval == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
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

######################################################################
## code for Mehta hypermutator data.
## modification of calc.cumulative.muts; that function specifically
## works on the LTEE metagenomics data.
## To normalize, we need to supply the number of sites at risk
## (such as sum of gene length).
## If plot.to.end is TRUE, then add one final row.
calc.Mehta.cumulative.muts <- function(d, d.metadata, plot.to.end=TRUE) {

    cumsum.per.pop.helper.func <- function(pop) {

        final.day <- 28 ## no more than 28 days in this experiment.
        
        df <- d %>% filter(Population==pop)
        if (nrow(df) == 0) { ## if no mutations in this pop.
            almost.done.df <- tibble(Population = pop,
                                     Day = final.day,
                                     count=0,
                                     cs=0)
        } else {
            summary.df <- df %>%
                arrange(t0) %>%
                group_by(Population,Day) %>%
                summarize(count=n()) %>%
                mutate(cs=cumsum(count)) %>%
                ungroup()
            ## if the final generation is not in ret.df,
            ## then add one final row (for nicer plots).
            final.row.df <- tibble(Population = pop,
                                   Day = final.day,
                                   count=max(summary.df$count),
                                   cs=max(summary.df$cs))
            if (plot.to.end) {
                almost.done.df <- bind_rows(summary.df, final.row.df)
            } else {
                almost.done.df <- summary.df
            }
        }
        ## add an row for Day == 0 (for nicer plots).
        init.row.df <- tibble(
            Population = pop,
            Day = 0,
            count = 0,
            cs = 0)
        
        ret.df <- bind_rows(init.row.df,almost.done.df)
        return(ret.df)
    }
    
    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my.genes <- d.metadata %>% dplyr::select(Gene, gene_length) %>% distinct()
    normalization.constant <- sum(my.genes$gene_length)
    
    c.dat <- map_dfr(.x=levels(d$Population),
                     .f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    
    return(c.dat)
}

## calc.traj.pvals, adapted for Mehta dataset.
calc.Mehta.traj.pvals <- function(data, PAO11.genes, gene.vec, N=10000, normalization.constant=NA) {

    ## check type for gene.vec (e.g., if factor, change to vanilla vector of strings)
    gene.vec <- as.character(gene.vec)
    
    ## each sample has the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## IMPORTANT: this function depends on variables defined in
    ## calculate.trajectory.tail.probs.
    generate.Mehta.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(PAO11.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(PAO11.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.Mehta.cumulative.muts(mut.subset, mut.subset.metadata,
                                                   normalization.constant) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }
        
    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),
                                         .f=generate.Mehta.cumulative.mut.subset)
    
    gene.vec.data <- data %>% filter(Gene %in% gene.vec)
    data.trajectory <- calc.Mehta.cumulative.muts(gene.vec.data,normalization.constant)
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
        group_by(Population,.drop=FALSE) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
    
    return(uppertail.probs)
}


########### Plotting code for Mehta hypermutator data.
## main difference is using Day instead of Generation,
## and calling calc.Mehta.cumulative.muts().
plot.Mehta.base.layer <- function(data, PAO11.genes, subset.size=50, N=1000, alphaval = 0.05, normalization.constant=NA, my.color="gray") {
    
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## IMPORTANT: this function depends on variables defined in plot.base.layer.
    generate.Mehta.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(PAO11.genes$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        mut.subset.metadata <- filter(PAO11.genes, Gene %in% rando.genes)
        c.mut.subset <- calc.Mehta.cumulative.muts(mut.subset, mut.subset.metadata,
                                                   normalization.constant) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),
                                         .f=generate.Mehta.cumulative.mut.subset)

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test. default is alphaval == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)

    p <- ggplot(filtered.trajectories,aes(x=Day,y=normalized.cs)) +
        ylab('Cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Day') +
        xlim(0, 28) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14)) +
        scale_y_continuous(labels=fancy_scientific,
                           breaks = scales::extended_breaks(n = 6),
                           limits = c(0, NA))
    return(p)
}


add.Mehta.cumulative.mut.layer <- function(p, layer.df, my.color) {
    p <- p +
        geom_point(data=layer.df,
                   aes(x=Day,y=normalized.cs),
                   color=my.color, size=0.2) +
        geom_step(data=layer.df, aes(x=Day,y=normalized.cs),
                  size=0.2, color=my.color)
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

## This plot visualizes a two-tailed test (alphaval = 0.05)
## against a bootstrapped null distribution for the derivative of cumulative mutations.
## This is what we want to use for publication.
plot.slope.of.base.layer <- function(data, subset.size=300, N=1000, alphaval = 0.05, normalization.constant=NA) {

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations,
    ## and then its derivative.
    generate.slope.of.cumulative.mut.subset <- function(idx) {
        rando.genes <- base::sample(unique(data$Gene),subset.size)
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

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories from each population,
    ## for a two-sided test. default is alphaval == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.D.norm.cs=max(D.normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.D.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
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
