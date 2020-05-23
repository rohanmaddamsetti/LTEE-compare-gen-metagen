## mutation-rate-analysis.R by Rohan Maddamsetti

## load library functions.
source("metagenomics-library.R")
library(viridis)
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
    ## This is for getting the colors right.
    mutate(Annotation=factor(Annotation)) %>% 
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
## for this reason, ALWAYS use mutation.data when examining mutation bias
## over the genome.

## gene.mutation.data only has mutations that have gene-level annotations.
## when examining genes on different strands, use gene.mutation.data.
############################

## we need a consistent color scale for all 6 classes of mutations in all plots.
## let's use the viridis color scheme.
## https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin

mut.annotation.vec <- levels(mutation.data$Annotation)
pal <- viridisLite::viridis(length(mut.annotation.vec))
names(pal) <- mut.annotation.vec
COL_SCALE <- scale_fill_manual(name = "Annotation", values = pal)
###################################################################################
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

make.summed.plot <- function(df, number.of.bins = 46) {
    ggplot(df, aes(x=Mbp.coordinate, fill=Annotation)) +
        geom_histogram(bins = number.of.bins) + 
        theme_classic() +
        ylab("Count") +
        xlab("Genomic position (Mb)") +
        theme(legend.position="bottom")
}

make.facet.mut.plot <- function(df) {
    make.summed.plot(df) + facet_wrap(.~Population,scales="free",nrow=4) 
}

point.mut.data <- mutation.data %>%
    filter(Annotation %in% c('missense', 'noncoding', 'nonsense', 'synonymous'))
indel.mut.data <- mutation.data %>% filter(Annotation %in% c('indel'))
sv.mut.data <- mutation.data %>% filter(Annotation %in% c('sv'))

## Figure 1: point mutations over the genome.
point.mut.plot <- make.facet.mut.plot(point.mut.data) + COL_SCALE
ggsave("../results/mutation-bias/figures/Fig1.pdf",point.mut.plot,width=7,height=7)

## S1 Figure: indels over the genome.
indel.mut.plot <- make.facet.mut.plot(indel.mut.data) +
    COL_SCALE +
    guides(fill=FALSE)
ggsave("../results/mutation-bias/figures/S1Fig.pdf",indel.mut.plot,width=6,height=4)

## S2 Figure: sv over the genome
sv.mut.plot <- make.facet.mut.plot(sv.mut.data) +
    COL_SCALE +
    guides(fill=FALSE)
ggsave("../results/mutation-bias/figures/S2Fig.pdf",sv.mut.plot,width=6,height=4)

## S3 Figure:
## all mutations from all pops summed over the genome.
## turn off legend so that x-axis is uniform across panels.
summed.point.mut.plot <- make.summed.plot(point.mut.data) +
    COL_SCALE +
    guides(fill=FALSE)
summed.indel.mut.plot <- make.summed.plot(indel.mut.data) +
    COL_SCALE +
    guides(fill=FALSE)
summed.sv.mut.plot <- make.summed.plot(sv.mut.data) +
    COL_SCALE +
    guides(fill=FALSE) 
summed.plot <- plot_grid(summed.point.mut.plot,
                         summed.indel.mut.plot,
                         summed.sv.mut.plot,
                         labels=c('A','B','C'),
                         ncol=1)
ggsave("../results/mutation-bias/figures/S3Fig.pdf",summed.plot,width=6,height=8)

## plot gene length against location in oriC coordinates.
gene.length.location.plot <- ggplot(REL606.genes, aes(x=oriC_start,y=gene_length)) +
    geom_point(size=0.5) + geom_smooth() + theme_classic()
ggsave("../results/mutation-bias/figures/gene_length_vs_location.pdf", gene.length.location.plot)
cor.test(REL606.genes$oriC_start, REL606.genes$gene_length, method="kendall")

###################################################################################
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

## 2) estimating local mutation rates by splitting the genome into z bins.
araplus3.local.probability.DNAtopology.not.hit <- partial(
    local.probability.that.DNAtopology.not.hit,
    araplus3.mut.data)

## 46 bins gives ~100 kB bins.
bins.to.try <- c(1,23,46,92,115,230,460,920)
map(bins.to.try,araplus3.local.probability.DNAtopology.not.hit)

#####################################################################################
## EVIDENCE OF STRAND-SPECIFIC BIAS.

## There is evidence of a strand-specific mutation bias on genes
## in the LTEE. Do genes on the lagging strand have a different number of
## mutations compared to the lagging strand?
## Will have to normalize by number of genes in each class.

## The key is that lagging/leading strands flip at the replication origin.
## so the asymmetry on each strand over the origin SHOWS the strand-specfic bias.
## the ratio of total mutations per strand on each side of the origin should give
## an estimate of the strength of this bias.

S4Fig <- make.summed.plot(gene.mutation.data) +
    facet_grid(strand~.,scales="fixed") +
    COL_SCALE
ggsave("../results/mutation-bias/figures/S4Fig.pdf", S4Fig, width=6, height=4)

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

########################################################################################
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

indel.plot <- plot.cumulative.muts(c.indel, my.color = pal[['indel']]) +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4)

sv.plot <- plot.cumulative.muts(c.sv, my.color = pal[['sv']]) +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4)

## using golden ratio phi for height/width ratio
Fig3 <- plot_grid(point.mut.plot,indel.plot, sv.plot,labels=c('A','B','C'),nrow=1,rel_heights=2/(1+sqrt(5)))
ggsave(filename="../results/mutation-bias/figures/Fig3.pdf",Fig3,width=7)
########################################################################################
## examine DNA repair and DNA polymerase/replication genes for mutator and anti-mutator
## candidates.

interesting.loci <- read.csv("../data/DNA-repair-and-replication.csv",header=TRUE,as.is=TRUE)

interesting.muts <- filter(mutation.data, Gene %in% unique(interesting.loci$Gene)) %>%
    filter(Annotation %in% c('missense','sv','indel','nonsense'))

interesting.parallel.nuc <- interesting.muts %>%
    group_by(Gene, Position) %>% summarize(count=n()) %>%
    arrange(desc(count)) %>% filter(count>1)

## Look for mutations in LTEE similar to those in:
## Deatherage et al. (2018). These are candidate anti-mutator alleles.

PResERV.genes <- c('polA','rne')
preserv.muts <- filter(mutation.data,Gene %in% PResERV.genes)
preserv.parallel <- preserv.muts %>% group_by(Gene,Position) %>%
    summarize(count=n()) %>% arrange(desc(count))

## There is one parallel mutation in rne: Position 1158096.
test <- filter(mutation.data,Position==1158096)

########################################################################################
## reanalysis of synonymous variation in natural populations.
## looks like problems with the K-S test... 

## goal was to investigate dS (and other classes of mutations) across the genome
## in the metagenomics data, revamping code from my 2015 Mol. Biol. Evol. paper.

## IMPORTANT ISSUE: changing the CDF (by changing rank on x-axis),
## while keeping the same set of probability masses for each gene
## changes K-S test statistics! How do I interpret this?
## Perhaps because K-S is for continuous and not discrete distributions?

## Right now I don't know the answer to this question. Report this as
## a Technical Comment to my 2015 MBE paper,
## in the Supplementary Text and in Supplementary Figure S1.

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
## set use.maddamsetti parameter to use one or other set of parameter estimates.
## for experimenting with how the ranking on the x-axis changes the statistics,
## I added the rev parameter to reverse the order of genes on the x-axis.
thetaS.KS.analysis <- function(the.data, REL606.genes, rank_by="length",use.maddamsetti=TRUE) {

    ## rank_by can have three values: "length", "thetaS", or "oriC".
    stopifnot(rank_by %in% c("length","thetaS","oriC"))
    
    ## 1) make an empirical cdf of mutations per core gene.
    ## do K-S tests for goodness of fit of the empirical cdf with cdfs for
    ## thetaS.
    
    hit.genes.df <- the.data %>%
        group_by(locus_tag, Gene, gene_length, oriC_start) %>%
        summarize(hits=n()) %>%
        ungroup()
    ## have to do it this complicated way, so that zeros are included.
    hit.genes.df <- full_join(REL606.genes,hit.genes.df)
    if (use.maddamsetti) {
        hit.genes.df <- mutate(hit.genes.df, thetaS=Maddamsetti_thetaS)
    } else {
        hit.genes.df <- mutate(hit.genes.df, thetaS=Martincorena_thetaS)
    }
    hit.genes.df <- hit.genes.df %>%
    filter(!(is.na(thetaS))) %>% ## only keep core genes.
        replace_na(list(hits=0)) %>%
        mutate(density=hits/gene_length)

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
                                  density=hit.genes.df$density,
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
        theme_classic()
    
    if (rank_by == "oriC") {
        ## find the index for gidA and mioC, which sandwich oriC.
        ## use these to plot the location of oriC.
        gidA.index <- filter(the.results.to.plot,Gene=='gidA')$index
        mioC.index <- filter(the.results.to.plot,Gene=='mioC')$index
        plot <- plot + scale_x_continuous('Genes ranked by chromosomal location') + geom_vline(xintercept=(gidA.index+mioC.index)/2,linetype='dotted',color='grey')
    } else if (rank_by == "length") {
        plot <- plot + scale_x_continuous('Genes ranked by length')
    } else if (rank_by == "thetaS") {
        plot <- plot + scale_x_continuous(expression(Genes~ranked~by~theta[s]))
        
    } else {
        stop("error 2 in make.KS.Figure.")
    }
    
    return(plot)
}

####################################################
gene.dS.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='synonymous')

gene.dN.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='missense')

## NOTE: thetaS KS analysis will only work for core genes!
## IMPORTANT!! RESULTS DEPEND ON X-AXIS ORDERING!
## This will be reported as a technical comment on my 2015 MBE paper.

cumsum.dS.core1 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"thetaS")
cumsum.dS.core2 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"length")
cumsum.dS.core3 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"oriC")

cumsum.dN.core1 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"thetaS")
cumsum.dN.core2 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"length")
cumsum.dN.core3 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"oriC")

dS.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dS.core1,"thetaS") + ggtitle("Synonymous mutations")
dS.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dS.core2,"length") + ggtitle("Synonymous mutations") 
dS.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dS.core3,"oriC") + ggtitle("Synonymous mutations") 

dN.thetaS.plot1 <- make.thetaS.KS.Figure(cumsum.dN.core1,"thetaS") + ggtitle("Missense mutations")
dN.thetaS.plot2 <- make.thetaS.KS.Figure(cumsum.dN.core2,"length") + ggtitle("Missense mutations")
dN.thetaS.plot3 <- make.thetaS.KS.Figure(cumsum.dN.core3,"oriC") + ggtitle("Missense mutations")

FigS10 <- plot_grid(dS.thetaS.plot1, dS.thetaS.plot2, dS.thetaS.plot3,
                   dN.thetaS.plot1, dN.thetaS.plot2, dN.thetaS.plot3,
                   labels=c('A','B','C','D','E','F'), nrow=2)
ggsave("../results/mutation-bias/figures/FigS10.pdf", height=7, width=10)

## Use pcor.R for partial correlation tests, freely available from:
## http://www.yilab.gatech.edu/pcor.R
source("pcor.R")

## Now, look at correlations between thetaS and dS.density and dN.density.
## marginally insignificant correlation using my estimates.
cor.test(cumsum.dS.core1$density,cumsum.dS.core1$thetaS, method="kendall")
## gene length still correlates with dS
cor.test(cumsum.dS.core1$gene_length,cumsum.dS.core1$density, method="kendall")
## gene length also correlates with thetaS
cor.test(cumsum.dS.core1$gene_length,cumsum.dS.core1$thetaS, method="kendall")

cor.test(cumsum.dN.core1$density,cumsum.dN.core1$thetaS, method="kendall")

cumsum.dS.core4 <- thetaS.KS.analysis(gene.dS.mutation.data,REL606.genes,"thetaS",use.maddamsetti=FALSE)
cumsum.dN.core4 <- thetaS.KS.analysis(gene.dN.mutation.data,REL606.genes,"thetaS",use.maddamsetti=FALSE)

## not significant at all when using Martincorena estimates.
cor.test(cumsum.dS.core4$density,cumsum.dS.core4$thetaS, method="kendall")
## even though gene length still correlates with dS and thetaS.
cor.test(cumsum.dS.core4$gene_length,cumsum.dS.core4$thetaS, method="kendall")
cor.test(cumsum.dS.core4$gene_length,cumsum.dS.core4$density, method="kendall")
