## metagenomics-analysis.R by Rohan Maddamsetti.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## TODO: do hierarchical modeling of point mutations in genes based on:
## 1) gene length
## 2) distance from peak in Ara+3.
## Question is to quantify contribution of mutation rate variation.

## TODO: Add Sastry modules to the analysis.

## TODO: measure the concentration/density of mutation as an empirical CDF.

## TODO: MAKE SURE MASKED REGIONS ARE TREATED PROPERLY!

## TODO: It could be quite interesting to examine genes in prophage that
## Jeff Barrick has annotated as misc_regions.

library(tidyverse)

##########################################################################
## FUNCTIONS
##########################################################################
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
        my.genes <- data %>% dplyr::select(Gene,gene_length) %>% distinct()
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
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops))) %>%
    ## This is for plotting mutation biases around oriC.
    mutate(oriC.coordinate=rotate.REL606.chr(Position,"oriC")) %>%
    mutate(terB.coordinate=rotate.REL606.chr(Position,"terB"))

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

## Now filter for nonsense, sv, and indels in genes (of course nonsense are always in genes)
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

## we treat nonsense as a nonsynonymous mutation-- for comparison to Tenaillon 2016 data.
gene.nodNdS.mutation.data <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    filter(Annotation %in% c("sv", "indel", "noncoding"))
##########################################################################################
## make base plot of null distributions by subsampling/bootstrapping.
## I will overlay data on top of these plots.

## Base plots here of null distributions: add the data lines on top to compare.
all.rando.plot <- plot.base.layer(gene.mutation.data, logscale=FALSE)
nodNdS.rando.plot <- plot.base.layer(gene.nodNdS.mutation.data, logscale=FALSE)
sv.indel.nonsen.rando.plot <- plot.base.layer(sv.indel.nonsense.gene.mutation.data, logscale=FALSE)

##########################################################################################
## MUTATION BIAS ANALYSIS.

## for indels and structural variation, we cannot distinguish between
## mutation hotspots vs. selection.

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

## For supplement, also plot with terminus as the origin in genome coordinates.
terB.mut.plot <- ggplot(mutation.data,aes(x=terB.coordinate,fill=Annotation)) +
    geom_histogram(bins=100) + 
    theme_classic() +
    facet_wrap(.~Population,scales="free") 
ggsave("../results/figures/terB-mutation-bias-histogram.pdf",terB.mut.plot,width=11,height=8)

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

##########################################################################################
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

## Figure 1 will be a conceptual figure of the different evolutionary models
## for different sets of genes.
## (purifying selection, coupon collecting, neutral or mutation accumulation, etc.)

## Figure 2 TODO: 
## combine these results to show distribution of different classes of mutations
## over the genome.

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
    dplyr::select(Gene,locus_tag,blattner,gene_length,product,nonsense.indel.sv.mut.density,nonsense.indel.sv.mut.count)


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
top_dN_checkme <- top.by.dN.density.genes %>% dplyr::select(Gene,locus_tag,blattner, gene_length,
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
    dplyr::select(-JW_id)

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
#########################################################################
## look at the accumulation of stars over time for top genes in the
## Tenaillon et al. (2016) genomics data.

## These curves are so far above the average, that comparing slopes is
## challenging. To solve this issue let's compare rates of mutation
## by looking at derivative itself.

## Look at these pictures carefully. Some indicate positive selection,
## while others indicate mutation accumulation.

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

## also look at the derivative of the accumulation of stars over time.
## base plot of null distribution for comparison.
D.of.all.rando.plot <- plot.slope.of.base.layer(gene.mutation.data)
D.of.nodNdS.rando.plot <- plot.slope.of.base.layer(gene.nodNdS.mutation.data)


## 1) plot top genes in non-mutators.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv')
top.nonmut.genomics <- top_n(nonmut.genomics, 50, wt=G.score)

top.nonmut.mutation.data <- gene.mutation.data %>%
                                   filter(Gene %in% top.nonmut.genomics$Gene.name)
c.top.nonmut.muts <- calc.cumulative.muts(top.nonmut.mutation.data)
D.of.c.top.nonmut.muts <- calc.slope.of.cumulative.muts(c.top.nonmut.muts)

tenaillon.plot1 <- all.rando.plot %>%
    add.cumulative.mut.layer(c.top.nonmut.muts,my.color="black",logscale=FALSE)

## evidence of continual fine tuning!
D.of.tenaillon.plot1 <- D.of.all.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.nonmut.muts,
                                      my.color='red')

## 2a) plot top genes in non-mutators, after excluding dN and dS mutations.
nonmut.nodNdS <- read.csv('../data/tenaillon2016-nonmutator-parallelism-nodSdN.csv')
top.nonmut.nodNdS <- top_n(nonmut.nodNdS, 50) ## using excluding dN dS column by default.

top.nonmut.nodNdS.mutation.data <- gene.mutation.data %>%
                                   filter(Gene %in% top.nonmut.nodNdS$Gene.name)
c.top.nonmut.nodNdS.muts <- calc.cumulative.muts(top.nonmut.nodNdS.mutation.data)
D.of.c.top.nonmut.nodNdS.muts <- calc.slope.of.cumulative.muts(c.top.nonmut.nodNdS.muts)

tenaillon.plot2a <- all.rando.plot %>%
    add.cumulative.mut.layer(c.top.nonmut.nodNdS.muts,my.color="black",logscale=FALSE)

D.of.tenaillon.plot2a <- D.of.all.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.nonmut.nodNdS.muts,
                                      my.color='red')


## 2b) exclude dN and dS in random comparison.
top.nonmut.nodNdS.mutation.data.2b <- gene.nodNdS.mutation.data %>%
                                   filter(Gene %in% top.nonmut.nodNdS$Gene.name)
c.top.nonmut.nodNdS.muts.2b <- calc.cumulative.muts(top.nonmut.nodNdS.mutation.data.2b)

tenaillon.plot2b <- nodNdS.rando.plot %>%
    add.cumulative.mut.layer(c.top.nonmut.nodNdS.muts.2b,my.color="black",logscale=FALSE)

D.of.tenaillon.plot2b <- D.of.nodNdS.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.nonmut.nodNdS.muts,
                                      my.color='red')


## 3) plot top genes in hypermutators.
mut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv')
top.mut.genomics <- top_n(mut.genomics, 50, wt=G.score)

top.mut.mutation.data <- gene.mutation.data %>%
    filter(Gene %in% top.mut.genomics$Gene.name)

c.top.mut.muts <- calc.cumulative.muts(top.mut.mutation.data)
D.of.c.top.mut.muts <- calc.slope.of.cumulative.muts(c.top.mut.muts)

tenaillon.plot3 <- all.rando.plot %>%
    add.cumulative.mut.layer(c.top.mut.muts,my.color="black",logscale=FALSE)

D.of.tenaillon.plot3 <- D.of.all.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.mut.muts,
                                      my.color='red')


## 4a) plot top genes in hypermutators, after excluding dN and dS mutations.

mut.nodNdS <- read.csv('../data/tenaillon2016-mutator-parallelism-nodSdN.csv')
top.mut.nodNdS <- top_n(mut.nodNdS, 50) ## using excluding dN dS column by default.

top.mut.nodNdS.mutation.data <- gene.mutation.data %>%
                                   filter(Gene %in% top.mut.nodNdS$Gene.name)
c.top.mut.nodNdS.muts <- calc.cumulative.muts(top.mut.nodNdS.mutation.data)
D.of.c.top.mut.nodNdS.muts <- calc.slope.of.cumulative.muts(c.top.mut.nodNdS.muts)

tenaillon.plot4a <- all.rando.plot %>%
    add.cumulative.mut.layer(c.top.mut.nodNdS.muts,my.color="black",logscale=FALSE)

D.of.tenaillon.plot4a <- D.of.all.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.mut.nodNdS.muts,my.color="red")


## 4b) exclude dN and dS in random comparison.
top.mut.nodNdS.mutation.data.4b <- gene.nodNdS.mutation.data %>%
                                   filter(Gene %in% top.mut.nodNdS$Gene.name)
c.top.mut.nodNdS.muts.4b <- calc.cumulative.muts(top.mut.nodNdS.mutation.data.4b)

tenaillon.plot4b <- nodNdS.rando.plot %>%
    add.cumulative.mut.layer(c.top.mut.nodNdS.muts.4b,my.color="black",logscale=FALSE)

D.of.tenaillon.plot4b <- D.of.nodNdS.rando.plot %>%
    add.slope.of.cumulative.mut.layer(D.of.c.top.mut.nodNdS.muts,my.color="red")



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

## plot for proteome sectors.
sector.plot <- all.rando.plot %>%
    add.cumulative.mut.layer(c.A.muts, my.color="black", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange", logscale=FALSE)
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

## now plot sv, indel, nonsense mutations on top.
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

eigen.plot <- all.rando.plot %>%
    add.cumulative.mut.layer(c.eigen1.muts, my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black", logscale=FALSE)
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

## why does Ara+3 show apparent selection in eigengene 6?
## Seems to be entirely driven by an excess of mutations in entF.
## at a first glance, maybe caused by gene conversion or homologous recombination
## with some other locus in the genome?
## If so, all these mutations should be linked in a cohort in the metagenomics data.
## This is not the case! Looks like a bona fide example of historical contingency!
eigengene6.mut.data %>% filter(Population=='Ara+3') %>%
    group_by(Gene) %>% summarize(count=n())

## plot eigen sv, indels, nonsense.
eigen.sv.indel.nonsen.plot <- sv.indel.nonsen.rando.plot %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen1.muts,my.color="red", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen2.muts,my.color="orange", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen3.muts,my.color="yellow", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen4.muts,my.color="green", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen5.muts,my.color="cyan", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen6.muts,my.color="blue", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen7.muts,my.color="violet", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen8.muts,my.color="pink", logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsen.eigen9.muts,my.color="black", logscale=FALSE)
ggsave(eigen.sv.indel.nonsen.plot,filename='../results/figures/eigen-sv-indel-nonsen-plot.png')

##########################################################################
## TODO:

## look at accumulation of stars over time for genes in different transcriptional
## modules inferred by Sastry et al. (2020) paper from Bernhard Palsson's group.


##########################################################################
## overall, no evidence of selection on dS.
## randomization tests show neutrality cannot be rejected when
## looking at dS density per gene (p = 0.14).

## to what extent is there evidence of intragenomic recombination between partially replicated chromosomes? intragenomic recombination/gene conversion?
## candidates for gene conversion: multiple dS in the same Population, in the same gene,
## and in the same cohort of mutations.
dS.in.same.cohort <- gene.dS.mutation.data %>% group_by(Population,Gene,t0,tf) %>% summarize(cohort=n()) %>% arrange(desc(cohort)) %>% filter(cohort>1) %>%
    filter(Gene!='ydfQ')

## might these be larger recombination events? Examine all kinds of mutations, where
## multiple mutations in the same gene in the same cohort.
multiple.muts.in.same.cohort <- gene.mutation.data %>% group_by(Population,Gene,t0,tf) %>% summarize(cohort=n()) %>% filter(cohort>1) %>% arrange(Population,t0,tf) %>% arrange(Gene)

## parallelism in dS at the same position in the same population
## we get exactly one gene: ydfQ--
## BUT THIS IS A BUG TO BE FIXED, CAUSED BY TWO LOCI WITH THE SAME BLATTNER NUMBER.
bug.to.fix <- gene.dS.mutation.data %>% group_by(Population,Gene,Position) %>% summarize(count=n()) %>% arrange(desc(count)) %>% filter(count>1)
buggy.ydfQ.mutations <- gene.dS.mutation.data %>% filter(Gene=='ydfQ')

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
