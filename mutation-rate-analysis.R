## mutation-rate-analysis.R by Rohan Maddamsetti

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## NOTE: Perhaps should cite my STLE paper-- recombination events tend to happen flanking
## the replication origin. Is that result connected to the wave-like mutation bias
## pattern seen here and in Patricia Foster's and Vaughn Cooper's evolution experiments?

## Also check out Hi-C papers and others reporting 3D chromosome in E. coli.

## load library functions.
source("metagenomics-library.R")
library(RColorBrewer)

##########################################################################
## MUTATION RATE AND BIASES DATA ANALYSIS
##########################################################################

## import thetaS estimates.
thetaS.estimates <- read.csv("../data/Maddamsetti2015_thetaS_estimates.csv")

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
    mutate(terB.coordinate=rotate.REL606.chr(Position,"terB")) %>%
    mutate(Mbp.coordinate=oriC.coordinate/1000000)

gene.mutation.data <- inner_join(mutation.data,REL606.genes)
## It turns out that some gene names map to multiple genes!!!
## There are only 4 such cases (8 genes). So just omit these from the analysis.
duplicate.genes <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    group_by(Gene) %>%
    summarize(checkme=length(unique(gene_length))) %>%
    filter(checkme>1)
## filter those duplicates.
gene.mutation.data <- gene.mutation.data %>%
    filter(!(Gene %in% duplicate.genes$Gene))

############################
### VERY IMPORTANT NOTE:
## mutation.data includes all mutations in the LTEE metagenomics data.
## for this reason, ALWAYS use mutation.data when examining mutation bias over the genome.

## gene.mutation.data only has mutations that have gene-level annotations.
## when examining genes on different strands, use gene.mutation.data.

##########################################################################################
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

make.summed.plot <- function(df) {
    ggplot(df, aes(x=Mbp.coordinate, fill=Annotation)) +
        ## use RColorBrewer scale
        scale_fill_brewer(palette = "RdYlBu", direction=-1,drop=FALSE) +
        geom_histogram(bins=100) + 
        theme_classic() +
        ylab("Count") +
        xlab("Genomic position (Mb)") +
        theme(legend.position="bottom")
}

make.facet.mut.plot <- function(df) {
    make.summed.plot(df) + facet_wrap(.~Population,scales="free") 
}

point.mut.data <- mutation.data %>%
    filter(Annotation %in% c('missense', 'synonymous', 'nonsense', 'noncoding'))
indel.mut.data <- mutation.data %>% filter(Annotation %in% c('indel'))
sv.mut.data <- mutation.data %>% filter(Annotation %in% c('sv'))

## Figure 1: point mutations over the genome.
point.mut.plot <- make.facet.mut.plot(point.mut.data)
ggsave("../results/mutation-bias/figures/Fig1.pdf",point.mut.plot,width=7,height=7)


## S1 Figure: indels over the genome.
indel.mut.plot <- make.facet.mut.plot(indel.mut.data) +
    scale_fill_manual(values="purple") + guides(fill=FALSE)
ggsave("../results/mutation-bias/figures/S1Fig.pdf",indel.mut.plot,width=6,height=4)

## S2 Figure: sv over the genome
sv.mut.plot <- make.facet.mut.plot(sv.mut.data) +
    scale_fill_manual(values="green") + guides(fill=FALSE)
ggsave("../results/mutation-bias/figures/S2Fig.pdf",sv.mut.plot,width=6,height=4)

## S3 Figure:
## all mutations from all pops summed over the genome.
## turn off legend so that x-axis is uniform across panels.
summed.point.mut.plot <- make.summed.plot(point.mut.data) + guides(fill=FALSE)
summed.indel.mut.plot <- make.summed.plot(indel.mut.data) +
    guides(fill=FALSE) +
    scale_fill_manual(values="purple")
summed.sv.mut.plot <- make.summed.plot(sv.mut.data) +
    guides(fill=FALSE) +
    scale_fill_manual(values="green")
summed.plot <- plot_grid(summed.point.mut.plot,
                         summed.indel.mut.plot,
                         summed.sv.mut.plot,
                         labels=c('A','B','C'),
                         ncol=1)
ggsave("../results/mutation-bias/figures/S3Fig.pdf",summed.plot,width=6,height=8)


## plot gene length against location in oriC coordinates.
gene.length.location.plot <- ggplot(REL606.genes, aes(x=oriC_start,y=gene_length)) +
    geom_point(size=0.5) + geom_smooth() + theme_classic()
gene.length.location.plot

##########################################################################################
### EPISTASIS AND HISTORICAL CONTINGENCY IN DNA TOPOLOGY GENES topA, fis, dusB.

## Calculate approximate probability of
## no mutations in fis, topA, dusB (yhdG) in Ara+3.

araplus3.mut.data <- filter(mutation.data,Population=='Ara+3')
topA.info <- filter(REL606.genes, Gene=='topA')
fis.info <- filter(REL606.genes, Gene=='fis')
dusB.info <- filter(REL606.genes, Gene=='yhdG')

## 1) assuming a uniform mutation rate.
araplus3.n.muts <- nrow(araplus3.mut.data)
fis.topA.dusB.length <- sum(c(topA.info$gene_length,
                              fis.info$gene_length,
                              dusB.info$gene_length))
uniform.probability.that.not.hit(fis.topA.dusB.length,araplus3.n.muts)

## 2) estimating local mutation rates by splitting the genome into z bins.

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

local.probability.that.DNAtopology.not.hit <- function(mutdata, z) {
    ## mutdata is ara+3 mutation data.
    ## z is the number of bins over the genome.
    ## THIS CODE ONLY WORKS FOR topA, fis, dusB.

    GENOME.LENGTH <- 4629812
    c <- GENOME.LENGTH/z ## length of each bin
    
    mutations.per.bin <- bin.mutations(araplus3.mut.data,z)
    dusB.bin <- find.bin(dusB.info, z)
    fis.bin <- find.bin(fis.info, z)
    topA.bin <- find.bin(topA.info, z)
    
    target.lengths <- c(dusB.info$gene_length,
                        fis.info$gene_length,
                        topA.info$gene_length)
    
    target.bin.mut.counts <- c(mutations.per.bin[dusB.bin],
                               mutations.per.bin[fis.bin],
                               mutations.per.bin[topA.bin])
    
    target.muts.per.base <- target.bin.mut.counts/c

    ## this is a dot product.
    ## multiply the muts.per.base by gene length, and sum.
    expected.num.muts.for.target <- sum(target.muts.per.base*target.lengths)

    total.muts <- nrow(araplus3.mut.data)
    expected.p.hit <- expected.num.muts.for.target/total.muts

    p.not.hit <- (1 - expected.p.hit)^total.muts
    return(p.not.hit)
}

araplus3.local.probability.DNAtopology.not.hit <- partial(
    local.probability.that.DNAtopology.not.hit,
    araplus3.mut.data)

## 46 bins gives ~100 kB bins.
bins.to.try <- c(1,23,46,92,115,230,460,920)
map(bins.to.try,araplus3.local.probability.DNAtopology.not.hit)

##########################################################################################
## EVIDENCE OF STRAND-SPECIFIC BIAS.

## There is evidence of a strand-specific mutation bias on genes
## in the LTEE. Do genes on the lagging strand have a different number of
## mutations compared to the lagging strand?
## Will have to normalize by number of genes in each class.

## The key is that lagging/leading strands flip at the replication origin.
## so the asymmetry on each strand over the origin SHOWS the strand-specfic bias.
## the ratio of total mutations per strand on each side of the origin should give
## an estimate of the strength of this bias.

FigS4 <- make.summed.plot(gene.mutation.data) +
    facet_grid(strand~.,scales="fixed") +
    scale_fill_brewer(palette = "RdPu",direction=-1,drop=FALSE)
ggsave("../results/mutation-bias/figures/FigS4.pdf", FigS4,width=6,height=4)

## Variable names assume that the lagging strand has more mutations--
## I can't tell based on arbitrary orientation labeling conventions for
## oriC.coordinate and strand! This doesn't matter in terms of reporting the results,
## however, as we measure asymmetry rather than the difference in number of
## mutations.

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

##########################################################################################
## MUTATION RATE ANALYSIS

##  plot of the cumulative number of different classes of mutations
## across the whole genome, in order to show evolution of mutation rates in LTEE.

## THESE RESULTS ARE REALLY COOL! WE SEE DIFFERENT KINDS OF ANTI-MUTATOR ALLELES!

## set normalization constant to 1, since we're looking over the whole genome.
point.mutation.data <- mutation.data %>%
    filter(Annotation %in% c("missense", "synonymous", "noncoding", "nonsense"))
sv.mutation.data <- mutation.data %>%
    filter(Annotation=='sv')
indel.mutation.data <- mutation.data %>%
    filter(Annotation=='indel')

## IMPORTANT: set plot.to.end parameter to FALSE.
c.point.muts <- calc.cumulative.muts(point.mutation.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)

c.sv <- calc.cumulative.muts(sv.mutation.data,
                             normalization.constant=1,
                             plot.to.end=FALSE)

c.indel <- calc.cumulative.muts(indel.mutation.data,                                       
                                normalization.constant=1,
                                plot.to.end=FALSE)

point.mut.plot <- plot.cumulative.muts(c.point.muts) +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4)

indel.plot <- plot.cumulative.muts(c.indel,my.color="purple") +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4)

sv.plot <- plot.cumulative.muts(c.sv,my.color="green") +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4)

## using golden ratio phi for height/width ratio
Fig3 <- plot_grid(point.mut.plot,indel.plot, sv.plot,labels=c('A','B','C'),nrow=1,rel_heights=2/(1+sqrt(5)))
ggsave(filename="../results/mutation-bias/figures/Fig3.pdf",Fig3,width=7)

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

    ## we don't plot alt2 for clarity (it doesn't change the message).
    plot <- ggplot(the.results.to.plot, aes(x=index)) +
        geom_line(aes(y=empirical), colour="red") + 
        geom_line(aes(y=null), linetype=2) +
        geom_line(aes(y=alt1), linetype='dotted') +
        scale_y_continuous('Cumulative proportion of mutations',limits=c(0,1)) +
        theme_classic() #+
        #theme(axis.title=element_text(size=18),axis.text=element_text(size=12))

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
## This will be reported as a technical comment on my 2015 MBE paper.

cumsum.dS.core1 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"thetaS")
cumsum.dS.core2 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"length")
cumsum.dS.core3 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"oriC")

cumsum.dN.core1 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"thetaS")
cumsum.dN.core2 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"length")
cumsum.dN.core3 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"oriC")

dS.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dS.core1,"thetaS")
dS.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dS.core2,"length")
dS.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dS.core3,"oriC")

dN.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dN.core1,"thetaS")
dN.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dN.core2,"length")
dN.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dN.core3,"oriC")

FigS5 <- plot_grid(dS.thetaS.plot1, dS.thetaS.plot2, dS.thetaS.plot3,
                   dN.thetaS.plot1, dN.thetaS.plot2, dN.thetaS.plot3,
                   labels=c('A','B','C','D','E','F'), nrow=2)
ggsave("../results/mutation-bias/figures/FigS5.pdf", height=7, width=10)


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
