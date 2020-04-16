## gene-module-analysis.R by Rohan Maddamsetti

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

source("metagenomics-library.R")

##########################################################################
## GENE MODULE DATA ANALYSIS
##########################################################################

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
    mutate(oriC_end=rotate.REL606.chr(end,"oriC"))

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

## There are four such cases: alr, bioD, ddl, maf (2 each)
## filter them from the analysis.
gene.mutation.data <- gene.mutation.data %>%
    filter(!(Gene %in% duplicate.genes$Gene))

##########################################################################
## BUT THIS IS A BUG TO BE FIXED, CAUSED BY TWO LOCI WITH THE SAME BLATTNER NUMBER.
bug.to.fix <- gene.mutation.data %>% group_by(Population,Gene,Position) %>% summarize(count=n()) %>% arrange(desc(count)) %>% filter(count>1)

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

## TODO: CHECK FOR BUGS IN MERGE. CHECK ydfQ.

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
## look at accumulation of stars over time for genes in different transcriptional
## modules inferred by Sastry et al. (2020) paper from Bernhard Palsson's group.

## Questions:
## Q1) are the regulators of the I-modulons under stronger selection than
## the genes that they regulate?

## Two ways of testing:
## A) group all regulators together, against all regulated genes.
## B) paired comparison. compare regulator against regulated genes for each I-modulon.

## Q2) Do any I-modulons show evidence of positive or purifying selection?
## I may have to do some kind of FDR correction for multiple hypothesis testing.

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

######################################################################################
