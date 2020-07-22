## mutation-rate-analysis.R by Rohan Maddamsetti

## load library functions.
source("metagenomics-library.R")
library(viridis)

##########################################################################
## FUNCTIONS
##########################################################################
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

## helper to map Allele to class of point mutation.
SNPToSpectrumMap <- function(SNP) {
    if (SNP == 'A->G') {
        return("A:T→G:C")
    } else if (SNP == "T->C") {
        return("A:T→G:C")
    } else if (SNP == "G->A") {
        return("G:C→A:T")
    } else if (SNP == "C->T") {
        return("G:C→A:T")
    } else if (SNP == "A->T") {
        return("A:T→T:A")
    } else if (SNP == "T->A") {
        return("A:T→T:A")
    } else if (SNP == "G->T") {
        return("G:C→T:A")
    } else if (SNP == "C->A") {
        return("G:C→T:A")
    } else if (SNP == "A->C") {
        return("A:T→C:G")
    } else if (SNP == "T->G") {
        return("A:T→C:G")
    } else if (SNP == "G->C") {
        return("G:C→C:G")
    } else if (SNP == "C->G") {
        return("G:C→C:G")
    }
}

SNPSpectrumToClassMap <- function(Spec) {
    if (Spec == 'A:T→G:C') {
        return("Transition")
    } else if (Spec == 'G:C→A:T') {
        return("Transition")
    } else {
        return("Transversion")
    }
}

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
    mutate(gene_length = strtoi(gene_length)) %>%
    mutate(oriC_start = rotate.REL606.chr(start,"oriC")) %>%
    mutate(oriC_end = rotate.REL606.chr(end,"oriC")) %>%
    ## annotate gene orientation for strand-specific bias analysis.
    mutate(gene_orientation = ifelse(
    ((oriC_start > 0) & (strand == 1)) | ((oriC_start < 0) & (strand == -1)),
    'L','R')) %>%
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

## in some analyses, we just want to look at hypermutator data.
hypermut.mutation.data <- mutation.data %>%
    filter(Population %in% hypermutator.pops)

## in some analyses, we just want to look at mutations in genes.
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

## in some analyses, we want to separately examine
## mismatch repair versus oxidative damage mutators.
## See Couce et al. (2017) in PNAS for this annotation.
MMR.mutator.pops <- c("Ara+3", "Ara-4", "Ara-3","Ara-2")
mutT.mutator.pops <- c("Ara-1", "Ara+6")

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

################################################################################
## MUTATION RATE ANALYSIS

##  plot of the cumulative number of different classes of mutations
## across the whole genome, in order to show evolution of mutation rates in LTEE.

## THESE RESULTS ARE REALLY COOL! WE SEE DIFFERENT KINDS OF ANTI-MUTATOR ALLELES!

## set normalization constant to 1, since we're looking over the whole genome.
point.mutation.data <- mutation.data %>%
    filter(Annotation %in% c("missense", "synonymous", "noncoding", "nonsense")) %>%
    mutate(Spectrum = sapply(Allele, SNPToSpectrumMap))

indel.mutation.data <- mutation.data %>%
    filter(Annotation=='indel')
sv.mutation.data <- mutation.data %>%
    filter(Annotation=='sv')
IS.sv.data <- sv.mutation.data %>%
    filter(str_detect(Allele, 'IS'))
sv.noIS.data <- sv.mutation.data %>%
    filter(str_detect(Allele, 'IS',negate = TRUE))

## IMPORTANT: set plot.to.end parameter to FALSE.
c.point.muts <- calc.cumulative.muts(point.mutation.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)

c.indel <- calc.cumulative.muts(indel.mutation.data,                                       
                                normalization.constant=1,
                                plot.to.end=FALSE)

c.IS.sv <- calc.cumulative.muts(IS.sv.data,
                             normalization.constant=1,
                             plot.to.end=FALSE)
c.sv.noIS <- calc.cumulative.muts(sv.noIS.data,
                                  normalization.constant=1,                                  
                                  plot.to.end=FALSE) 

point.mut.plot <- plot.cumulative.muts(c.point.muts,my.color = 'red') +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4) +
    ggtitle("Point mutations")

indel.plot <- plot.cumulative.muts(c.indel, my.color = pal[['indel']]) +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4) +
    ggtitle("Indel mutations")

sv.plot <- plot.cumulative.muts(c.IS.sv, my.color = pal[['sv']]) %>%
    add.cumulative.mut.layer(c.sv.noIS, my.color = "gray") +
    ylab("Cumulative number of mutations") +
    facet_wrap(.~Population,scales='fixed',nrow=4) +
    ggtitle("Structural variants")

## using golden ratio phi for height/width ratio
Fig1 <- plot_grid(point.mut.plot,indel.plot, sv.plot,labels=c('A','B','C'),nrow=1,rel_heights=2/(1+sqrt(5)))
ggsave(filename="../results/mutation-bias/figures/Fig1.pdf",Fig1,width=7)

################################################################################
## Figures 2 and 3: time course trajectory plots for
## oxidative damage repair and mismatch repair genes.
################################################################################
## EVIDENCE OF STRAND-SPECIFIC BIAS.

## Null comparison: plot nucleotide composition on the sense strand
## for each gene, based on its gene orientation.

nuc.composition.df <- REL606.genes %>%
    gather(`G`,`C`,`A`,`T`,key="Nucleotide",value="NucleotideCount") %>%
    mutate(oriC.coordinate = (oriC_start + oriC_end)/2) %>%
    mutate(Mbp.coordinate=oriC.coordinate/1000000) %>%
    ## hack to make the plot I want:
    ## make a big df with a row for every nucleotide.
    uncount(NucleotideCount)

Fig4A <- ggplot(nuc.composition.df, aes(x=Mbp.coordinate, fill=Nucleotide)) +
    geom_histogram(bins = 46) + ## TODO: change to parameter 'number.of.bins'
    theme_classic() +
    ggtitle("Null expectation based on\nnucleotide composition") +
    ylab("Count") +
    xlab("Genomic position (Mb)") +
    theme(legend.position="bottom") +
    facet_grid(strand~., scales="fixed") +
    scale_fill_viridis_d(option = "inferno") +
    geom_vline(xintercept=0,color='grey',linetype='dotted')

## There is evidence of a strand-specific mutation bias on genes
## in the LTEE. Do genes on the lagging strand have a different number of
## mutations compared to the lagging strand?
## Will have to normalize by number of genes in each class.

## The key is that lagging/leading strands flip at the replication origin.
## so the asymmetry on each strand over the origin SHOWS the strand-specfic bias.
## the ratio of total mutations per strand on each side of the origin should give
## an estimate of the strength of this bias.

## The statistics indicate strand-specific mutation biases for both MMR and MutT
## hypermutators-- however those biases go in opposite directions!
## By examining asymmetries in nucleotide composition (i.e. GC skew),
## it seems that asymmetric context-dependent biases are probably the
## best explanation (Way Sung et al. 2015 in Molecular Biology and Evolution).

hypermut.dS.data <- gene.mutation.data %>%
    filter(Annotation == "synonymous") %>%
    filter(Population %in% hypermutator.pops) %>%
    mutate(MMR.deficient = ifelse(Population %in% MMR.mutator.pops,
                                  "MMR hypermutators",
                                  "MutT hypermutators")) %>%
    mutate(Spectrum = sapply(Allele, SNPToSpectrumMap)) %>%
    mutate(SNPClass = sapply(Spectrum, SNPSpectrumToClassMap))

Fig4B <- ggplot(hypermut.dS.data,
               aes(x=Mbp.coordinate, fill=Spectrum)) +
    geom_histogram(bins = 46) + ## TODO: change to parameter 'number.of.bins'
    theme_classic() +
    ylab("Count") +
    xlab("Genomic position (Mb)") +
    theme(legend.position="bottom") +
    ggtitle("Observed mutation distribution\nper strand") +
    facet_grid(strand~MMR.deficient, scales="fixed") +
    scale_fill_viridis_d(option = "plasma") +
    geom_vline(xintercept=0,color='grey',linetype='dotted')

Fig4 <- plot_grid(Fig4A, Fig4B, labels= c('A','B'),nrow=1)
Fig4 ## I save this image using quartz() because it renders the arrows properly.
## I am not sure how to set up the proper device backend for cowplot/ggsave
## to render the arrows properly.

## I call one gene orientation 'L' and the other gene orientation 'R'.
## The fact that these labels are arbitrary don't matter in terms of
## reporting the results because we measure asymmetry
##rather than the difference in number of mutations.

MMR.L.gene.mut.data <- hypermut.dS.data %>%
    filter(gene_orientation == 'L') %>%
    filter(MMR.deficient == "MMR hypermutators")

MMR.R.gene.mut.data <- hypermut.dS.data %>%
    filter(Population %in% hypermutator.pops) %>%
    filter(gene_orientation == 'R') %>%
    filter(MMR.deficient == "MMR hypermutators")

MutT.L.gene.mut.data <- hypermut.dS.data %>%
    filter(gene_orientation == 'L') %>%
    filter(MMR.deficient == "MutT hypermutators")

MutT.R.gene.mut.data <- hypermut.dS.data %>%
    filter(Population %in% hypermutator.pops) %>%
    filter(gene_orientation == 'R') %>%
    filter(MMR.deficient == "MutT hypermutators")

R.REL606.genes <- REL606.genes %>% filter(gene_orientation == 'R')
L.REL606.genes <- REL606.genes %>% filter(gene_orientation == 'L')
nrow(R.REL606.genes)
nrow(L.REL606.genes)

nrow(MMR.R.gene.mut.data)
nrow(MMR.L.gene.mut.data)

nrow(MutT.R.gene.mut.data)
nrow(MutT.L.gene.mut.data)

R.gene.target <- sum(R.REL606.genes$gene_length)
L.gene.target <- sum(L.REL606.genes$gene_length)

MMR.strand.mut.vec <- c(nrow(MMR.R.gene.mut.data),
                        nrow(MMR.L.gene.mut.data))
binom.test(x = MMR.strand.mut.vec, p = R.gene.target/(R.gene.target+L.gene.target))

MutT.strand.mut.vec <- c(nrow(MutT.R.gene.mut.data),
                         nrow(MutT.L.gene.mut.data))
binom.test(x=MutT.strand.mut.vec,p=R.gene.target/(R.gene.target+L.gene.target))

##################################################################################
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

## Figure 5: mutations over the genome in hypermutator populations.
hypermut.plot <- make.facet.mut.plot(hypermut.mutation.data) + COL_SCALE
ggsave("../results/mutation-bias/figures/Fig5.pdf", hypermut.plot,width=7,height=7)

no.Araplus3.data <- hypermut.mutation.data %>%
    filter(Population != "Ara+3") %>%
    mutate(MMR.deficient = ifelse(Population %in% MMR.mutator.pops,
                                  "MMR hypermutators (excluding Ara+3)",
                                  "MutT hypermutators"))

## Fig 6. When excluding Ara+3, the MMR mutators show a weaker wave.
## the mutT mutators, on the other hand, don't show the wave.
Fig6 <- make.summed.plot(no.Araplus3.data) + COL_SCALE +
    facet_wrap(.~MMR.deficient)

ggsave("../results/mutation-bias/figures/Fig6.pdf",Fig6,width=6,height=3.5)

###################################################################################
### EPISTASIS AND HISTORICAL CONTINGENCY IN DNA TOPOLOGY GENES topA, fis, dusB.

## Calculate approximate probability of
## no mutations in fis, topA, dusB (yhdG) in Ara+3.
## Ara+3 has one synonymous mutation in yhdG, that
## went to fixation in the cohort with the mutS/mutH mutator alleles,
## around 4000 generations.
## SO EXCLUDE SYNONYMOUS MUTATIONS FOR THIS CALCULATION!

## but maybe something interesting with synonymous mutations in dusB/yhdG?
## Figure S4 is intriguing... looks like synonymous mutations in dusB are
## overdispersed (i.e. one in each lineage).
dusB.muts <- mutation.data %>% filter(Gene == 'yhdG') %>%
    arrange(Position)
## YES! There is parallelism at the exact mutation in Ara+3! Also shows
## up early in Ara+6.

araplus3.mut.data <- filter(mutation.data,Population=='Ara+3') %>%
    filter(Annotation != "synonymous")
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

#################################
## is there purifying selection on synonymous mutations in topA, dusB, and fis?
dS.mut.data <- gene.mutation.data %>%
    filter(Annotation=='synonymous')
DNA.topology.dS.muts <- dS.mut.data %>%
    filter(Gene %in% c('topA','fis','yhdG'))

## run STIMS to test for purifying selection on dS in topA, fis, and yhdG.
## highly significant signal of purifying selection on dS in Ara-1 and Ara-3.
DNA.topology.dS.pvals <- calc.traj.pvals(dS.mut.data, c('topA','fis','yhdG'))

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

#################################################################################
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

araplus3.m3 <- glm(hits ~ 0 + thetaS, family="poisson", data=araplus3.hit.genes.df)
summary(araplus3.m3) ## worse than araplus3.m1

#################################################################################
## TODO: Use this Supplementary Figure as a main figure. !!!!!!!!!!!!!!!!!
## Supplementary Figure 1: break down Figure 1A by mutation spectrum.

hypermut.SNPs <- mutation.data %>%
    filter(Population %in% hypermutator.pops) %>%
    ## relevel for plotting.
    mutate(Population = factor(Population)) %>%
    filter(Annotation %in% c("synonymous","missense","nonsense","noncoding")) %>%
    mutate(Spectrum = sapply(Allele, SNPToSpectrumMap))

SNPclass1.data <- hypermut.SNPs %>% filter(Spectrum == "A:T→C:G")
SNPclass2.data <-  hypermut.SNPs %>% filter(Spectrum == "A:T→G:C")
SNPclass3.data <- hypermut.SNPs %>% filter(Spectrum == "A:T→T:A")
SNPclass4.data <- hypermut.SNPs %>% filter(Spectrum == "G:C→A:T")
SNPclass5.data <- hypermut.SNPs %>% filter(Spectrum == "G:C→C:G")
SNPclass6.data <- hypermut.SNPs %>% filter(Spectrum == "G:C→T:A")

SNP.annotation.vec <- sort(unique(hypermut.SNPs$Spectrum))
SNP.pal <- viridisLite::viridis(length(SNP.annotation.vec),option="plasma")
names(SNP.pal) <- SNP.annotation.vec
SNP_COL_SCALE <- scale_fill_manual(name = "Spectrum", values = SNP.pal)

c.SNP1 <- calc.cumulative.muts(SNPclass1.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)
c.SNP2 <- calc.cumulative.muts(SNPclass2.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)
c.SNP3 <- calc.cumulative.muts(SNPclass3.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)
c.SNP4 <- calc.cumulative.muts(SNPclass4.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)
c.SNP5 <- calc.cumulative.muts(SNPclass5.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)
c.SNP6 <- calc.cumulative.muts(SNPclass6.data,
                                     normalization.constant=1,
                                     plot.to.end=FALSE)

S1FigA <- plot.cumulative.muts(c.SNP1, my.color = SNP.pal[['A:T→C:G']]) %>%
    add.cumulative.mut.layer(c.SNP2, my.color = SNP.pal[['A:T→G:C']]) %>%
    add.cumulative.mut.layer(c.SNP3, my.color = SNP.pal[['A:T→T:A']]) %>%
    add.cumulative.mut.layer(c.SNP4, my.color = SNP.pal[['G:C→A:T']]) %>%
    add.cumulative.mut.layer(c.SNP5, my.color = SNP.pal[['G:C→C:G']]) %>%
    add.cumulative.mut.layer(c.SNP6, my.color = SNP.pal[['G:C→T:A']]) + 
    ylab("Cumulative number of mutations") +
    ggtitle("Spectrum of point mutations over time in hypermutator populations") +
    facet_wrap(.~Population,scales='free',nrow=1) +
    theme(axis.title.y = element_text(size = 10))

S1FigB <- plot.cumulative.muts(c.SNP1, my.color = SNP.pal[['A:T→C:G']]) %>%
    add.cumulative.mut.layer(c.SNP2, my.color = SNP.pal[['A:T→G:C']]) %>%
    add.cumulative.mut.layer(c.SNP3, my.color = SNP.pal[['A:T→T:A']]) %>%
    add.cumulative.mut.layer(c.SNP4, my.color = SNP.pal[['G:C→A:T']]) %>%
    add.cumulative.mut.layer(c.SNP5, my.color = SNP.pal[['G:C→C:G']]) %>%
    add.cumulative.mut.layer(c.SNP6, my.color = SNP.pal[['G:C→T:A']]) + 
    ylab("Cumulative number of mutations") +
    ggtitle("Inset showing non-dominant point mutation spectra over time") +
    facet_wrap(.~Population,scales='free',nrow=1) +
    ylim(0,750) +
    theme(axis.title.y = element_text(size = 10))

## grab and use the figure legend from Figure 4B.
S1_legend <- cowplot::get_legend(Fig4B)
grid::grid.newpage()
grid::grid.draw(S1_legend)

S1Fig <- plot_grid(plot_grid(S1FigA,S1FigB, labels=c('A','B'),ncol=1),
                   plot_grid(S1_legend),ncol=1,rel_heights=c(5,0.5))
## plot.
S1Fig

#################################################################################
## Table 1: putative mutator and anti-mutator alleles.

Table1.data <- gene.mutation.data %>%
    select(Population, Position, Gene, Allele, Annotation, t0, tf,
           transit_time, fixation, final_frequency, product)

Ara.minus.1.data <- Table1.data %>%
    filter(Population=='Ara-1') %>%
    filter(Gene %in% c('uvrC','mutT','mutY')) %>%
    arrange(t0)
