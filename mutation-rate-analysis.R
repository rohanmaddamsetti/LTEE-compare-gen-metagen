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
Fig1 <- plot_grid(point.mut.plot,indel.plot, sv.plot,labels=c('A','B','C'),nrow=1,rel_heights=2/(1+sqrt(5)))
ggsave(filename="../results/mutation-bias/figures/Fig1.pdf",Fig3,width=7)
########################################################################################
## Figures 2 and 3: time course trajectory plots for
## oxidative damage repair and mismatch repair genes.
########################################################################################
## EVIDENCE OF STRAND-SPECIFIC BIAS.

## There is evidence of a strand-specific mutation bias on genes
## in the LTEE. Do genes on the lagging strand have a different number of
## mutations compared to the lagging strand?
## Will have to normalize by number of genes in each class.

## The key is that lagging/leading strands flip at the replication origin.
## so the asymmetry on each strand over the origin SHOWS the strand-specfic bias.
## the ratio of total mutations per strand on each side of the origin should give
## an estimate of the strength of this bias.

Fig4 <- make.summed.plot(gene.mutation.data) +
    facet_grid(strand~.,scales="fixed") +
    COL_SCALE
ggsave("../results/mutation-bias/figures/Fig4.pdf", Fig4, width=6, height=4)

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

hypermut.mutation.data <- mutation.data %>%
    filter(Population %in% hypermutator.pops)
## Figure 5: mutations over the genome in hypermutator populations.
hypermut.plot <- make.facet.mut.plot(hypermut.mutation.data) + COL_SCALE
ggsave("../results/mutation-bias/figures/Fig5.pdf", hypermut.plot,width=7,height=7)

## Fig S4. Excluding Ara+3, don't see any wave.
no.Araplus3 <- hypermut.mutation.data %>%
    filter(Population != 'Ara+3')
no.Araplus3.plot <- make.summed.plot(no.Araplus3) + COL_SCALE
ggsave("../results/mutation-bias/figures/FigS4.pdf", no.Araplus3.plot,width=4,height=4)

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
## examine DNA repair and DNA polymerase/replication genes for mutator and anti-mutator
## candidates.

interesting.loci <- read.csv("../data/DNA-repair-and-replication.csv",header=TRUE,as.is=TRUE)

interesting.muts <- filter(mutation.data, Gene %in% unique(interesting.loci$Gene)) %>%
    filter(Annotation %in% c('missense','sv','indel','nonsense'))

interesting.parallel.nuc <- interesting.muts %>%
    group_by(Gene, Position) %>% summarize(count=n()) %>%
    arrange(desc(count)) %>% filter(count>1)

most.interesting.muts <- filter(mutation.data,Position %in% interesting.parallel.nuc$Position)

########################################################################################
## reanalysis of synonymous variation in natural populations.
## looks like I set up the K-S test incorrectly in my
## my 2015 Mol. Biol. Evol. paper.

## this reanalysis finds the same conclusion. So it seems like my
## 2015 Mol. Bio. Evol. paper came to the "right" conclusion
## using dodgy statistics.

gene.dS.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='synonymous')

## analysis over all LTEE populations.
dS.summary.df <- gene.dS.mutation.data %>%
    group_by(locus_tag, Gene, gene_length, oriC_start) %>%
    summarize(hits=n()) %>%
    ungroup()
## have to do it this complicated way, so that zeros are included.
hit.genes.df <- full_join(REL606.genes,dS.summary.df) %>%
    ## use the corrected thetaS estimates in Martincorena et al. (2012).
    mutate(thetaS=Martincorena_thetaS) %>%
    ## treat thetaS as if its an estimate of per bp mutation rate,
    ## as advocated by Martincorena et al. (2012).
    mutate(thetaS_by_length=thetaS*gene_length) %>%
    filter(!(is.na(thetaS))) %>% ## only keep core genes.
    replace_na(list(hits=0)) %>%
    mutate(density=hits/gene_length)

## analysis over all LTEE populations.
## let's use a Poisson GLM approach.
m1 <- glm(hits ~ 0 + gene_length, family="poisson", data=hit.genes.df)
summary(m1)

m2 <- glm(hits ~ 0 + thetaS_by_length, family="poisson", data=hit.genes.df)
summary(m2) ## worse than m1

m3 <- glm(hits ~ 0 + thetaS, family="poisson", data=hit.genes.df)
summary(m3) ## worse than m1

## calculate for individual LTEE populations.
pop.dS.summary.df <- gene.dS.mutation.data %>%
    ## key difference is group_by Population here too.
    group_by(Population,locus_tag, Gene, gene_length, oriC_start) %>%
    summarize(hits=n()) %>%
    ungroup()
pop.hit.genes.df <- full_join(REL606.genes, pop.dS.summary.df) %>%
    mutate(thetaS=Martincorena_thetaS) %>%
    mutate(thetaS_by_length=thetaS*gene_length) %>%
    filter(!(is.na(thetaS))) %>% ## only keep core genes.
    replace_na(list(hits=0)) %>%
    mutate(density=hits/gene_length)

araplus3.hit.genes.df <- pop.hit.genes.df %>%
    filter(Population == 'Ara+3')

## let's use a Poisson GLM approach.
araplus3.m1 <- glm(hits ~ 0 + gene_length, family="poisson", data=araplus3.hit.genes.df)
summary(araplus3.m1)

araplus3.m2 <- glm(hits ~ 0 + thetaS_by_length, family="poisson", data=araplus3.hit.genes.df)
summary(araplus3.m2) ## worse than araplus3.m1

araplus3.m3 <- glm(hits ~ 0 + thetaS, family="poisson", data=ara.plus3.hit.genes.df)
summary(araplus3.m3) ## worse than araplus3.m1
