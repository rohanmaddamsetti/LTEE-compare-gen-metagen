## metagenomics-library.R by Rohan Maddamsetti
## This file is a library of common functions used in both paper 2A and paper 2B.

library(tidyverse)
library(cowplot)

#########################################################################
## A nice function for memory management/checking usage.
## From: https://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session/11899573
# improved list of objects
.ls.objects <- function(pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

##########################################################################
## FUNCTIONS FOR DATA ANALYSIS
##########################################################################

fancy_scientific <- function(x) {
    ## function for plotting better y-axis labels.
    ## see solution here for nice scientific notation on axes.
    ## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

uniform.probability.that.not.hit <- function(l,n) {
    ## calculate the probability that a locus (or set of loci) is not hit by mutations,
    ## assuming uniform mutation rate.
    ## l is locus length, and n is the number of mutations in the dataset.

    GENOME.LENGTH <- 4629812
    p <- (1 - (l/GENOME.LENGTH))^n
    return(p)
}

rotate.REL606.chr <- function(my.position, c) {
    #' function to rotate REL606 genome coordinates,
    #' setting oriC/terB at the center of plots
    #' that examine mutation bias over the chromosome.
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

#######################################

.construct.cumsum.per.pop.helper <- function(d, final.time, reset.pop.levels, pop) {
    ## This function is never used directly. It generates functions that
    ## calculate the cumulative number of mutations in LTEE data or Mehta
    ## data, depending on the setting of d, final.time, and reset.pop.levels.
    
    ## This is basically for the cases where
    ## only one population is being examined, or we don't want to reorder.
    pop.levels <- levels(d$Population)
    
    if (reset.pop.levels) { ## This is only useful for LTEE data.
        ## This constant is to make sure that all pops are in the levels
        ## of the Population factor after mergers, etc.
        pop.levels <- c("Ara-5","Ara-6", "Ara+1", "Ara+2",
                        "Ara+4", "Ara+5", "Ara-1", "Ara-2",
                        "Ara-3", "Ara-4", "Ara+3", "Ara+6")
    }
    
    df <- d %>% filter(Population==pop)
    
    if (nrow(df) == 0) { ## if no mutations in this pop.
        almost.done.df <- tibble(Population = factor(pop, levels = pop.levels),
                                 Time = final.time,
                                 count=0,
                                 cs=0)
    } else {
        summary.df <- df %>%
            group_by(Population, Time) %>%
            summarize(count=n(), .groups = "drop_last") %>%
            mutate(cs=cumsum(count)) %>%
            ungroup()
        ## if the final generation is not in ret.df,
        ## then add one final row (for nicer plots).
        final.row.df <- tibble(Population = factor(pop, levels = pop.levels),
                               Time = final.time,
                               count=0,
                               cs=max(summary.df$cs))   
        almost.done.df <- bind_rows(summary.df, final.row.df)
    }
    
    ## add an row for Time == 0 (for nicer plots).
    init.row.df <- tibble(
        Population = factor(pop, levels = pop.levels),
        Time = 0,
        count = 0,
        cs = 0)
    
    ret.df <- bind_rows(init.row.df,almost.done.df)
    return(ret.df)
}

.calc.cumulative.muts <- function(ltee.not.mehta, reset.pop.levels,
                                 d, d.metadata) {
    ## look at accumulation of stars over time
    ## in other words, look at the rates at which the mutations occur over time.
    ## To normalize, we need to supply the number of sites at risk
    ## (such as sum of gene length).
    ## If plot.to.end is TRUE, then add one final row.

    ## This function should not be run directly, rather, specialized versions of
    ## this are generated next, as syntactic sugar.
    
    if (ltee.not.mehta) {
        ## Generate a one-argument helper function for LTEE data.
        ## set the last timepoint to 6.3 * 10000 generations.
        cumsum.helper <- partial(
            .construct.cumsum.per.pop.helper,
            d, 6.3, reset.pop.levels)
    } else {
        ## Generate a one-argument helper function for Mehta data.
        ## set the last timepoint to 28 days and never reset levels(Population).
        cumsum.helper <- partial(
            .construct.cumsum.per.pop.helper,
            d, 28, FALSE)
    }

    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my.genes <- d.metadata %>% dplyr::select(Gene,gene_length) %>% distinct()
    normalization.constant <- sum(my.genes$gene_length)
    
    c.dat <- map_dfr(.x = levels(d$Population),
                     .f = cumsum.helper) %>%
        mutate(normalized.cs = cs / normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    
    return(c.dat)
}

## These two function are what are actually called in context.
calc.LTEE.cumulative.muts <- partial(.f = .calc.cumulative.muts,
                                     TRUE, TRUE)

calc.single.LTEE.pop.cumulative.muts <- partial(.f = .calc.cumulative.muts,
                                     TRUE, FALSE)

calc.Mehta.cumulative.muts <- partial(.f = .calc.cumulative.muts,
                                     FALSE, FALSE)
###############################
## the next two functions are for estimating local mutation rates by
## splitting the genome into z bins.

## IMPORTANT: THESE FUNCTIONS ONLY WORK ON REL606.

## find.bin is useful for sampling random genes while preserving bin identity,
## that is, only sampling genes near the genes of interest.

bin.mutations <- function(mut.data, z) {
    ## z is the number of bins we are using.
    ## count the number of mutations in each bin, M(z_i).
    ## filter mut.data using each adjacent pair of fenceposts,
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
        bin.data <- mut.data %>%
            filter(Position >= left & Position < right)
        bin.mut.count <- nrow(bin.data)
        mutations.by.bin.vec[i] <- bin.mut.count
    }
    
    ## assert that all mutations have been assigned to a bin.
    stopifnot(sum(mutations.by.bin.vec) == nrow(mut.data))
    
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

assign.genes.to.bins <- function(d.metadata) {
    ## use 46 bins so that each bin is ~10,000 bp (for REL606).
    find.bin.46 <- partial(find.bin,z=46)
    
    ## output a bin from each row in gene.info.
    bin.list <- d.metadata %>%
        split(.$Gene) %>%
        lapply(find.bin.46)
    
    ## associate each gene in the genome to its bin, and add as a column.
    bin.df <- data.frame(Gene = names(bin.list), bin = as.numeric(bin.list))
    return(bin.df)
}

sample.genes.by.genomebin <- function(module.metadata.with.bins) {
    ## map each gene in gene.vec to a random gene in its bin
    ## in the genome.
    sample_one <- partial(sample_n,size = 1)
    
    module.to.random.module <- module.metadata.with.bins %>%
        ## the genes to sample are stored in a nested column called "data".
        mutate(sampled.info = map(data, sample_one)) %>%
        transmute(Gene, sampled.gene = map_chr(sampled.info, function(x) x$Gene))
    return(module.to.random.module$sampled.gene)
}


calc.LTEE.cumulative.mut.subset.by.bin <- function(d, d.metadata,
                                                   module.metadata.with.bins,
                                                   idx) {
    ## This function is partialized within bootstrap.LTEE.traj.by.loc.
    rando.genes <- sample.genes.by.genomebin(module.metadata.with.bins)
    mut.subset <- filter(d, Gene %in% rando.genes)
    mut.subset.metadata <- filter(d.metadata, Gene %in% rando.genes)
    c.mut.subset <- calc.LTEE.cumulative.muts(mut.subset, mut.subset.metadata) %>%
        mutate(bootstrap_replicate=idx)
    return(c.mut.subset)
}

bootstrap.LTEE.traj.by.loc <- function(d, d.metadata, gene.vec, N) {
    ## sample genes near the genes in the module of interest,
    ## i.e. gene.vec.
    
    ## for efficiency, pre-calculate gene bin assignments before bootstrapping.
    bin.df <- assign.genes.to.bins(d.metadata)
    
    ## map bins to the genes in those bins (excepting those in gene.vec).
    nested.genes.to.sample.by.bin <- left_join(d.metadata, bin.df) %>%
        ## when sampling, exclude the genes in the module of interest.
        filter(!(Gene %in% gene.vec)) %>%
        group_by(bin) %>%
        nest()
    
    ## add the genes to be sampled as a list column corresponding
    ## to the genes in gene.vec,
    ## such that each gene in the module maps to the other genes in its bin.
    ## then, one of those genes will be sampled to make the random module,
    ## while preserving bin location.
    module.metadata.with.bins <- d.metadata %>%
        filter(Gene %in% gene.vec) %>%
        left_join(bin.df) %>%
        left_join(nested.genes.to.sample.by.bin)

    ## partialized function to calculate cumulative mutations.
    generate.LTEE.cumulative.mut.subset.by.bin <- partial(
        .f = calc.LTEE.cumulative.mut.subset.by.bin,
        d, d.metadata, module.metadata.with.bins)
    
    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped.trajectories <- map_dfr(
        .x=seq_len(N),
        .f=generate.LTEE.cumulative.mut.subset.by.bin)

    return(bootstrapped.trajectories)
}

filter.pop.trajectories.by.data <- function(data.trajectory, pop.trajectories) {

    data.trajectory.summary <- data.trajectory %>%
        group_by(Population, .drop=FALSE) %>%
        summarize(final.norm.cs = max(normalized.cs)) %>%
        ungroup()
    
    pop <- unique(pop.trajectories$Population)
    data.traj <- filter(data.trajectory.summary, Population == pop)
    
    final.data.norm.cs <- unique(data.traj$final.norm.cs)
    tail.trajectories <- filter(pop.trajectories, final.norm.cs >= final.data.norm.cs)
    return(tail.trajectories)
}

calc.traj.pvals <- function(d, d.metadata, gene.vec, N=10000,
                            ltee.not.mehta = TRUE, sample.genes.by.location=FALSE) {
    ## calculate the tail probabilities of the true cumulative mutation trajectory
    ## of a given vector of genes (a 'module'), based on resampling
    ## random sets of genes. Returns the upper tail of null distribution,
    ## or P(random trajectory >= the actual trajectory).
    ## Output: a dataframe with three columns: Population, count, p.val

    ## check type for gene.vec (e.g., if factor, change to vanilla vector of strings)
    gene.vec <- as.character(gene.vec)
    ## each sample has the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)

    gene.vec.data <- d %>% filter(Gene %in% gene.vec)
    gene.vec.metadata <- d.metadata %>% filter(Gene %in% gene.vec)

    if(ltee.not.mehta) {
        data.trajectory <- calc.LTEE.cumulative.muts(gene.vec.data, gene.vec.metadata)
    } else {
        data.trajectory <- calc.Mehta.cumulative.muts(gene.vec.data, gene.vec.metadata)
    }

    ## partialized-- the last blocks uses this function to calculate the tail.
    trajectory.filter.helper <- partial(
        .f = filter.pop.trajectories.by.data,
        data.trajectory)

    if(!ltee.not.mehta) { ## if analyzing the Mehta dataset:
        ## This function takes the index for the current draw, and samples the data,
        ## generating a random gene set for which to calculate cumulative mutations.
        generate.Mehta.cumulative.mut.subset <- partial(
            .f = .generate.cumulative.mut.subset,
            d, d.metadata, subset.size, calc.Mehta.cumulative.muts)
        
        ## make a dataframe of bootstrapped trajectories.
        ## look at accumulation of stars over time for random subsets of genes.
        bootstrapped.trajectories <- map_dfr(
            .x=seq_len(N),
            .f=generate.Mehta.cumulative.mut.subset)
        
    } else { ## then analyzing LTEE dataset.
        
        if (sample.genes.by.location) {
            ## then sample genes near the genes in the module of interest,
            ## i.e. gene.vec.
            bootstrapped.trajectories <- bootstrap.LTEE.traj.by.loc(
                d, d.metadata, gene.vec, N)
            
        } else {
            ## This function takes the index for the current draw, and samples the data,
            ## generating a random gene set for which to calculate cumulative mutations.
            generate.LTEE.cumulative.mut.subset <- partial(
                .f = .generate.cumulative.mut.subset,
                d, d.metadata, subset.size, calc.LTEE.cumulative.muts)
            
            ## make a dataframe of bootstrapped trajectories.
            ## look at accumulation of stars over time for random subsets of genes.
            bootstrapped.trajectories <- map_dfr(
                .x = seq_len(N),
                .f = generate.LTEE.cumulative.mut.subset)          
        }
    }
    
    uppertail.probs <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() %>%
        ## split by Population, then filter for bootstraps > data trajectory.
        split(.$Population) %>%
        map_dfr(.f=trajectory.filter.helper) %>%
        group_by(Population,.drop=FALSE) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
    
    return(uppertail.probs)
}

get.middle.trajectories <- function(bootstrapped.trajectories, alphaval) {
    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## usual default is alphaval == 0.05.
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
    return(filtered.trajectories)
}

get.middle.trajs.over.all.pops <- function(bootstrapped.trajectories, alphaval) {
    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## usual default is alphaval == 0.05.
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup()
    
    top.trajectories <- trajectory.summary %>%
        slice_max(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        slice_min(prop=alphaval/2, order_by=final.norm.cs) %>%
        dplyr::select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        dplyr::select(-in.top,-in.bottom)
    return(filtered.trajectories)
}

.generate.cumulative.mut.subset <- function(data, gene.metadata, subset.size,
                                            calc.cumulative.muts.function, idx) {
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    ## This function is not used directly-- it is 
    rando.genes <- base::sample(unique(gene.metadata$Gene), subset.size)
    mut.subset <- filter(data, Gene %in% rando.genes)
    mut.subset.metadata <- filter(gene.metadata, Gene %in% rando.genes)
    c.mut.subset <- calc.cumulative.muts.function(mut.subset, mut.subset.metadata) %>%
        mutate(bootstrap_replicate=idx)
    return(c.mut.subset)
}

plot.base.layer <- function(data, gene.metadata, subset.size=50, N=1000,
                            alphaval = 0.05, my.color = "gray", ltee.not.mehta=TRUE) {
    ## This plot visualizes a two-tailed test (alphaval = 0.05)
    ## against a bootstrapped null distribution.
    ## Throughout, plots use the minimum subsample size to subsample
    ## the null distribution, to increase the variance in order to
    ## make a conservative comparison.

    if (ltee.not.mehta) { ## use function for LTEE data.
        generate.cumulative.mut.subset <- partial(
            .f = .generate.cumulative.mut.subset,
            data, gene.metadata, subset.size, calc.LTEE.cumulative.muts)
    } else { ## use function for Mehta data.
        generate.cumulative.mut.subset <- partial(
            .f = .generate.cumulative.mut.subset,
            data, gene.metadata, subset.size, calc.Mehta.cumulative.muts)        
    }
    
    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped.trajectories <- map_dfr(
        .x=seq_len(N),
        .f=generate.cumulative.mut.subset)

    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## default is alphaval == 0.05.
    middle.trajs <- get.middle.trajectories(bootstrapped.trajectories, alphaval)
    
    p <- ggplot(middle.trajs,aes(x=Time,y=normalized.cs)) +
           ylab('Cumulative number of mutations (normalized)') +
           theme_classic() +
           geom_point(size=0.2, color=my.color) +
           theme(axis.title.x = element_text(size=14),
                 axis.title.y = element_text(size=14),
                 axis.text.x  = element_text(size=14),
                 axis.text.y  = element_text(size=14)) +
           scale_y_continuous(labels=fancy_scientific,
                              breaks = scales::extended_breaks(n = 6),
                              limits = c(0, NA)) +
        facet_wrap(.~Population, scales='free', nrow=4)
    
    if (ltee.not.mehta) {
        p <- p +
            xlab('Generations (x 10,000)') +
            xlim(0,6.3)
    } else {
        p <- p +
            xlab('Day') +
            xlim(0, 28)        
    }
    
    return(p)
}

add.LTEE.base.layer <- function(p, data, REL606.genes, my.color, subset.size=50,
                           N=1000, alphaval = 0.05) {
    ## add a base layer to a plot. used in Imodulon code.
    generate.cumulative.mut.subset <- partial(
        .f = .generate.cumulative.mut.subset,
        data, REL606.genes, subset.size, calc.LTEE.cumulative.muts)
    
    bootstrapped.trajectories <- map_dfr(
        .x=seq_len(N), .f=generate.cumulative.mut.subset)

    middle.trajs <- get.middle.trajectories(bootstrapped.trajectories, alphaval)

    p <- p + geom_point(data=middle.trajs,
                        aes(x=Time, y=normalized.cs),
                        size=0.2, color=my.color)
    return(p)                
}

add.cumulative.mut.layer <- function(p, layer.df, my.color) {
    ## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
    p <- p +
        geom_point(data = layer.df,
                   aes(x = Time, y = normalized.cs),
                   color = my.color, size = 0.2) +
        geom_step(data = layer.df, aes(x = Time, y = normalized.cs),
                  size = 0.2, color = my.color)
    return(p)
}

plot.cumulative.muts <- function(mut.data, my.color="black", ltee.not.mehta=TRUE) {
    ## calculate cumulative numbers of mutations in each category.
    ## for vanilla plotting, without null distributions, as plotted by
    ## plot.base.layer.
    p <- ggplot(mut.data, aes(x = Time, y = normalized.cs)) +
        ylab('Cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size = 0.2, color = my.color) +
        geom_step(size = 0.2, color = my.color) +
        facet_wrap(. ~ Population, scales = 'free', nrow = 4)
    
    if (ltee.not.mehta) {
        p <- p +
            xlab("Generations (x 10,000)") +
            xlim(0, 6.3)            
    } else {
        p <- p +
            xlab("Days") +
            xlim(0, 28)                    
    }
    
    return(p)
}

######################################################################
## versions of STIMS plotting code, summing over all populations.

calc.LTEE.cumulative.muts.over.all.pops <- function(d, d.metadata) {
    ## look at accumulation of stars over time
    ## in other words, look at the rates at which the mutations occur over time.
    ## To normalize, we need to supply the number of sites at risk
    ## (such as sum of gene length).
    ## If plot.to.end is TRUE, then add one final row.
    ## This function sums mutations over ALL populations.

    finalgen <- 6.3 ## this is outside of the data collection
    ## for nice plotting (final generation in mutation.data is 6.275).
    
    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my.genes <- d.metadata %>% dplyr::select(Gene,gene_length) %>% distinct()
    normalization.constant <- sum(my.genes$gene_length)
    
    summary.df <- d %>%
        group_by(Time) %>%
        summarize(count=n(), .groups = "drop_last") %>%
        mutate(cs=cumsum(count)) %>%
        ungroup()
    ## if the final generation is not in ret.df,
    ## then add one final row (for nicer plots).
    final.row.df <- tibble(Time=finalgen,
                           count=0,
                           cs=max(summary.df$cs))
    
    almost.done.df <- bind_rows(summary.df, final.row.df)

    ## add an row for Time == 0 (for nicer plots).
    init.row.df <- tibble(Time = 0,
                          count = 0,
                          cs = 0)
    
    c.dat <- bind_rows(init.row.df,almost.done.df) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    
    return(c.dat)
}

plot.base.layer.over.all.pops <- function(data, REL606.genes, subset.size=50,
                                          N=1000, alphaval = 0.05, my.color="gray") {
    ## This plot visualizes a two-tailed test (alphaval = 0.05)
    ## against a bootstrapped null distribution.
    ## Throughout, plots use the minimum subsample size to subsample
    ## the null distribution, to increase the variance in order to
    ## make a conservative comparison.

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset.over.all.pops <- partial(
        .f = .generate.cumulative.mut.subset,
        data, REL606.genes, subset.size, calc.LTEE.cumulative.muts.over.all.pops)

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.
    
    bootstrapped.trajectories <- map_dfr(
        .x = seq_len(N),
        .f = generate.cumulative.mut.subset.over.all.pops)

    middle.trajs <- get.middle.trajs.over.all.pops(
        bootstrapped.trajectories, alphaval)

    p <- ggplot(middle.trajs,aes(x=Time,y=normalized.cs)) +
        ylab('Cumulative number of mutations (normalized)') +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14)) +
        scale_y_continuous(labels=fancy_scientific,
                           breaks = scales::extended_breaks(n = 6),
                           limits = c(0, NA)) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3)
    
    return(p)
}

add.base.layer.over.all.pops <- function(p, data, REL606.genes, my.color, subset.size=50, N=1000, alphaval = 0.05) {
    ## add a base layer to a plot. used in Imodulon code.
    
    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset <- partial(
        .f = .generate.cumulative.mut.subset,
        data, REL606.genes, subset.size, calc.LTEE.cumulative.muts.over.all.pops)

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.
    
    bootstrapped.trajectories <- map_dfr(
        .x=seq_len(N),
        .f=generate.cumulative.mut.subset)

    middle.trajs <- get.middle.trajs.over.all.pops(bootstrapped.trajectories, alphaval)

    p <- p + geom_point(data = middle.trajs,
                        aes(x = Time, y = normalized.cs),
                        size = 0.2, color = my.color)
    return(p)                
}

calc.all.pops.LTEE.cumulative.mut.subset.by.bin <- function(d, d.metadata,
                                                            module.metadata.with.bins,
                                                            idx) {
    ## This function is partialized within bootstrap.all.LTEE.pops.traj.by.loc.
    rando.genes <- sample.genes.by.genomebin(module.metadata.with.bins)
    mut.subset <- filter(d, Gene %in% rando.genes)
    mut.subset.metadata <- filter(d.metadata, Gene %in% rando.genes)
    c.mut.subset <- calc.LTEE.cumulative.muts.over.all.pops(
        mut.subset, mut.subset.metadata) %>%
        mutate(bootstrap_replicate=idx)
    return(c.mut.subset)
}

bootstrap.all.LTEE.pops.traj.by.loc <- function(data, REL606.genes, gene.vec, N) {
    ## see bootstrap.LTEE.traj.by.loc. This version does not group by population.

    ## for efficiency, pre-calculate gene bin assignments before bootstrapping.
    bin.df <- assign.genes.to.bins(d.metadata)
    
    ## map bins to the genes in those bins (excepting those in gene.vec).
    nested.genes.to.sample.by.bin <- left_join(d.metadata, bin.df) %>%
        ## when sampling, exclude the genes in the module of interest.
        filter(!(Gene %in% gene.vec)) %>%
        group_by(bin) %>%
        nest()
    
    ## add the genes to be sampled as a list column corresponding
    ## to the genes in gene.vec,
    ## such that each gene in the module maps to the other genes in its bin.
    ## then, one of those genes will be sampled to make the random module,
    ## while preserving bin location.
    module.metadata.with.bins <- d.metadata %>%
        filter(Gene %in% gene.vec) %>%
        left_join(bin.df) %>%
        left_join(nested.genes.to.sample.by.bin)
    
    ## partialized function to calculate cumulative mutations.
    generate.all.LTEE.pops.cumulative.mut.subset.by.bin <- partial(
        .f = calc.all.pops.LTEE.cumulative.mut.subset.by.bin,
        d, d.metadata, module.info.with.bins)
    
    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    bootstrapped.trajectories <- map_dfr(
        .x=seq_len(N),
        .f=generate.all.LTEE.pops.cumulative.mut.subset.by.bin)

    return(bootstrapped.trajectories)
    
}

calc.all.pops.traj.pvals <- function(d, d.metadata, gene.vec, N=10000, sample.genes.by.location=FALSE) {
    ## calculate the tail probabilities of the true cumulative mutation trajectory
    ## of a given vector of genes (a 'module'), based on resampling
    ## random sets of genes. Returns the upper tail of null distribution,
    ## or P(random trajectory >= the actual trajectory).
    ## Output: a dataframe with two columns: count, p.val

    ## check type for gene.vec (e.g., if factor, change to vanilla vector of strings)
    gene.vec <- as.character(gene.vec)
    ## each sample has the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)

    gene.vec.data <- d %>% filter(Gene %in% gene.vec)
    gene.vec.metadata <- d.metadata %>% filter(Gene %in% gene.vec)

    data.trajectory <- calc.LTEE.cumulative.muts.over.all.pops(
        gene.vec.data, gene.vec.metadata)
    ## This is a number, not a vector.
    final.data.norm.cs <- max(data.trajectory$normalized.cs)


    if (sample.genes.by.location) {
        ## then sample genes near the genes in the module of interest,
        ## i.e. gene.vec.
        bootstrapped.trajectories <- bootstrap.LTEE.traj.by.loc(
            d, d.metadata, gene.vec, N)
        
        
        generate.all.LTEE.pops.cumulative.mut.subset.by.loc <- function(idx) {
            rando.genes <- sample.genes.by.genomebin()
            mut.subset <- filter(d, Gene %in% rando.genes)
            mut.subset.metadata <- filter(d.metadata, Gene %in% rando.genes)
            c.mut.subset <- calc.LTEE.cumulative.muts.over.all.pops(
                mut.subset,
                mut.subset.metadata) %>%
                mutate(bootstrap_replicate=idx)
            return(c.mut.subset)
        }
        
        ## make a dataframe of bootstrapped trajectories.
        ## look at accumulation of stars over time for random subsets of genes.
        bootstrapped.trajectories <- map_dfr(
            .x=seq_len(N),
            .f=generate.all.LTEE.pops.cumulative.mut.subset.by.loc)
        
    } else {
        ## This function takes the index for the current draw, and samples the data,
        ## generating a random gene set for which to calculate cumulative mutations.
        generate.all.pops.LTEE.cumulative.mut.subset <- partial(
            .f = .generate.cumulative.mut.subset,
            d, d.metadata, subset.size,
            calc.LTEE.cumulative.muts.over.all.pops)
        
        ## make a dataframe of bootstrapped trajectories.
        ## look at accumulation of stars over time for random subsets of genes.
        bootstrapped.trajectories <- map_dfr(
            .x = seq_len(N),
            .f = generate.all.pops.LTEE.cumulative.mut.subset)
    }
    
    uppertail.probs <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() %>%    
        ## filter for bootstraps > data trajectory.
        filter(final.norm.cs >= final.data.norm.cs) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
    
    return(uppertail.probs)
}

################################################################################
## Calculate densities of mutations per gene, over all LTEE and within each
## population.

calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    ## Examine the distribution of various classes of mutations across genes in the
    ## genomics or metagenomics data. Which genes are enriched? Which genes are depleted?
    ## Then, can look at the annotation of these genes in STRING.
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

pop.calc.gene.mutation.density <- function(gene.mutation.data, mut_type_vec) {
    ## Examine the gene mutation density by population.
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
## analysis of tempo of metagenomic evolution using mutations per unit time
## rather than cumulative number of mutations.
## maybe this is in some sense analogous to spike trains in neuroscience?
## maybe in the future techniques for analyzing neural spikes may be useful
## for analyzing spikes in observed mutations in a particular gene set over time?
## This code is not currently used for anything.

calc.slope.of.cumulative.muts <- function(c.muts) {
    ## Calculate the derivative of the cumulative accumulation of mutation occurrence.
    ## This is simply the rate of mutation occurrence in a class of genes.

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

plot.slope.of.cumulative.muts <- function(mut.data, my.color="black") {
    ## calculate derivative of cumulative numbers of mutations in each category.
    p <- ggplot(mut.data,aes(x=Time,y=D.normalized.cs)) +
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

add.slope.of.cumulative.mut.layer <- function(p, layer.df, my.color) {
    ## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
    p <- p +
        geom_point(data=layer.df, aes(x=Time,y=D.normalized.cs), color=my.color, size=0.2) +
        geom_step(data=layer.df, aes(x=Time,y=D.normalized.cs), color=my.color, size=0.2) +
        geom_smooth(data=layer.df,size=0.2, color=my.color)
    return(p)
}

plot.slope.of.base.layer <- function(data, subset.size=300, N=1000, alphaval = 0.05, normalization.constant=NA) {
    ## This plot visualizes a two-tailed test (alphaval = 0.05)
    ## against a bootstrapped null distribution for the derivative
    ## of cumulative mutations. This is what we want to use for publication.

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

    filtered.trajectories <- filter.trajectories(bootstrapped.trajectories, alphaval)

    p <- ggplot(filtered.trajectories,aes(x=Time,y=D.normalized.cs)) +
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
