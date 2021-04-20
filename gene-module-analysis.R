## gene-module-analysis.R by Rohan Maddamsetti

source("metagenomics-library.R")

##########################################################################
## GENE MODULE DATA ANALYSIS
##########################################################################

## get the lengths of all protein-coding genes in REL606.
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
## Some gene names map to multiple genes! Filter these out.
duplicate.genes <- REL606.genes %>%
    group_by(Gene) %>%
    summarize(copies=length(unique(gene_length))) %>%
    filter(copies>1)
## There are four such cases: alr, bioD, ddl, maf (2 each for 8 loci total).
## filter them from this analysis.
REL606.genes <- REL606.genes %>%
    filter(!(Gene %in% duplicate.genes$Gene))

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv
gene.mutation.data <- read.csv(
    '../results/LTEE-metagenome-mutations.csv',
    header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    ## This for changing the ordering of populations in plots.
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops))) %>%
    inner_join(REL606.genes) %>%
    filter(Gene!='intergenic')

## get mutation parallelism in the LTEE genomes published in Tenaillon et al. (2016).
## This is used for filtering essential genes
## in the purifying selection control analysis,
## and is used in the positive selection control analysis.
nonmut.genomics <- read.csv('../data/tenaillon2016-nonmutator-parallelism.csv') %>%
    ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)    

hypermut.genomics <- read.csv('../data/tenaillon2016-mutator-parallelism.csv') %>%
    ## make sure these genes passed the filters on REL606.genes.
    filter(Gene.name %in% REL606.genes$Gene)    

top.nonmut.genomics <- top_n(nonmut.genomics, 50, wt=G.score)
top.hypermut.genomics <- top_n(hypermut.genomics, 50, wt=G.score)

#########################################################################
## RELAXED SELECTION CONTROL EXPERIMENT.

## list of neutral genes from Alejandro Couce.
## His description: these are known 'neutral' genes
## we use to calibrate our analyses (i.e., genes we know to be cryptic, non
## expressed or have been experimentally shown to be neutral; see attached file).
## Note: this list of genes includes citT.

neutral.genes <- read.csv("../data/neutral_compilation.csv", header=TRUE,as.is=TRUE) %>%
    ## make sure that only loci that pass filters are included in the analysis.
    filter(Gene %in% REL606.genes$Gene)

neutral.genes.metadata <- REL606.genes %>% filter(Gene %in% neutral.genes$Gene)

neutral.mut.data <- gene.mutation.data %>%
    filter(Gene %in% neutral.genes$Gene)

c.neutral.genes <- calc.cumulative.muts(neutral.mut.data, neutral.genes.metadata)

neutral.base.layer <- plot.base.layer(
    gene.mutation.data,
    REL606.genes,
    subset.size=length(unique(neutral.genes$Gene)))

## Figure 2: plot of "gold standard" neutral genes.
Fig2 <- neutral.base.layer %>% 
    add.cumulative.mut.layer(c.neutral.genes, my.color="black")
ggsave("../results/gene-modules/figures/Fig2.pdf", Fig2)

## calculate more rigorous statistics than the figures.
neutral.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(neutral.genes$Gene))
## recalculate, sampling from the same genomic regions (bins).
## results should be unchanged.
neutral.pvals.loc <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(neutral.genes$Gene),sample.genes.by.location=TRUE)

#########################################################################
## PURIFYING SELECTION CONTROL EXPERIMENT AND ANALYSIS.

## Get essential and near-essential genes reported in
## Supplementary Table 1 of Couce et al. 2017.
## I manually fixed the names of a couple genes in this dataset.
## The original names are in the "Name" column, and updated names
## are in the "Gene" column.
essential.genes <- read.csv("../data/Couce2017-LTEE-essential.csv") %>%
    inner_join(REL606.genes) %>% filter(!(is.na(locus_tag)))

## a significant proportion of genes under positive selection in the LTEE are
## essential genes, as reported in Maddamsetti et al. (2017).
## filter these ones out, since we are interested in purifying selection.
## 21 out of 50 top non-mut genes are essential.
nonmut.top.hit.essential <- essential.genes %>%
    filter(Gene %in% top.nonmut.genomics$Gene.name)
## what about the hypermutators? 3 out of 50 top hypermut genes.
hypermut.top.hit.essential <- essential.genes %>%
    filter(Gene %in% top.hypermut.genomics$Gene.name)

## ## filtering out top G-score genes in the LTEE genomics dataset.
purifying.genes <- essential.genes %>%
    filter(!(Gene %in% nonmut.top.hit.essential$Gene)) %>%
    filter(!(Gene %in% hypermut.top.hit.essential$Gene))

purifying.mut.data <- gene.mutation.data %>%
    filter(Gene %in% purifying.genes$Gene)
c.purifying.genes <- calc.cumulative.muts(purifying.mut.data, purifying.genes)

purifying.base.layer <- plot.base.layer(
    gene.mutation.data,
    REL606.genes,
    subset.size=length(unique(purifying.genes$Gene)))

##  Yes. evidence of purifying selection in these genes based on my test.
Fig3 <- purifying.base.layer %>% 
    add.cumulative.mut.layer(c.purifying.genes, my.color="black")
ggsave("../results/gene-modules/figures/Fig3.pdf", Fig3)

## calculate more rigorous statistics than the figures.
purifying.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes, unique(purifying.genes$Gene))
## recalculate, sampling from the same genomic regions (bins).
## results should be unchanged.
purifying.pvals.loc <- calc.traj.pvals(gene.mutation.data, REL606.genes,
                                                       unique(purifying.genes$Gene),
                                                       sample.genes.by.location=TRUE)

## What is the association between genes that are essential by transposon
## mutagenesis, and genes that are completely depleted in mutations?
## Find all genes that have no mutations whatsoever in any population.
## IMPORTANT: most of these genes are very short. Many of them are probably
## depleted by chance, and are not significant after FDR-correction.
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

no.mutation.genes <- REL606.genes %>%
    ## no hits allowed in the metagenomics.
    filter(!(Gene %in% gene.mutation.data$Gene)) %>%
    ## and no hits (other than dS, amps, and intergenic) allowed in the genomics.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length)) %>%
    ## estimate the probability that these genes were not hit by chance.
    mutate(pval = uniform.probability.that.not.hit(gene_length, nrow(gene.mutation.data))) %>%
    mutate(fdr.qval = p.adjust(pval,"fdr"))

## Look at genes that only have dS.
all.mutation.density <- calc.gene.mutation.density(
    gene.mutation.data,
    c("missense", "sv", "synonymous", "noncoding", "indel", "nonsense")) %>%
    rename(all.mut.count = mut.count) %>%
    rename(all.mut.density = density)
all.except.dS.density <- calc.gene.mutation.density(
    gene.mutation.data,c("sv", "indel", "nonsense", "missense")) %>%
    rename(all.except.dS.mut.count = mut.count) %>%
    rename(all.except.dS.mut.density = density)
## combine these into one dataframe.
gene.mutation.densities <- REL606.genes %>%
    full_join(all.mutation.density) %>%
    full_join(all.except.dS.density)
#### CRITICAL STEP: replace NAs with zeros.
#### We need to keep track of genes that haven't been hit by any mutations
#### in a given mutation class (sv, indels, dN, etc.)
gene.mutation.densities[is.na(gene.mutation.densities)] <- 0
gene.mutation.densities <- tbl_df(gene.mutation.densities)

## need this next line to calculate p-values.
no.dS.count <- sum(all.except.dS.density$all.except.dS.mut.count)
## genes that are only affected by synonymous mutations.
only.dS.allowed.genes <- gene.mutation.densities %>%
    filter(all.except.dS.mut.count == 0) %>%
    ## no hits (other than dS, amps, and intergenic) allowed in the LTEE genomics
    ## to remove false positives.
    filter(!(str_detect(LTEE.genomics.mutated.genestr,Gene))) %>%
    arrange(desc(gene_length)) %>%
    ## estimate the probability that these genes were only hit by dS by chance.
    mutate(pval = uniform.probability.that.not.hit(gene_length,no.dS.count)) %>%
    mutate(fdr.qval = p.adjust(pval,"fdr"))

## 105 genes just have dS.
only.dS.genes <- only.dS.allowed.genes %>% filter(!(Gene %in% no.mutation.genes$Gene))
## 65 genes don't have even dS in the metagenomic data.
no.dS.genes <- only.dS.allowed.genes %>% filter(Gene %in% no.mutation.genes$Gene)

## If these sets are under purifying selection, then they should be more essential.
## let's examine overlap with Couce Tn10 mutagenesis report.
## make this into a table to report in the paper.

## 18 of these.
no.dS.genes.in.purifying.set <- no.dS.genes %>%
    filter(Gene %in% purifying.genes$Gene)
write.csv(no.dS.genes.in.purifying.set,file="../results/gene-modules/Table1.csv")

## 30 of these.
only.dS.genes.in.purifying.set <- only.dS.genes %>%
    filter(Gene %in% purifying.genes$Gene)
write.csv(only.dS.genes.in.purifying.set,file="../results/gene-modules/Table2.csv")

## the overlap is highly significant:
## 18 genes in purifying.genes and no.dS.genes.
## 64 no.dS genes.

## 30 genes in purifying.genes and only.dS.genes.
## 105 only.dS genes.

## 491 purifying.genes after filters.
## 3948 total genes after filters.

## contingency table for association between genes in Couce data and no dS in LTEE.
fisher.test(matrix(c(18,473,46,3411),2))

## contingency table for association between genes in Couce data and only dS in LTEE.
fisher.test(matrix(c(30,461,75,3382),2))

#################################################################################
## POSITIVE SELECTION CONTROL EXPERIMENT AND ANALYSIS.

## look at the accumulation of stars over time for top genes in the
## Tenaillon et al. (2016) genomics data.

## base plots of null distribution for comparison.
rando.plot <- plot.base.layer(gene.mutation.data, REL606.genes)
    
## 1) plot top genes in non-mutators.
top.nonmut.mutation.data <- gene.mutation.data %>%
    filter(Gene %in% top.nonmut.genomics$Gene.name)

c.top.nonmuts <- calc.cumulative.muts(top.nonmut.mutation.data, top.nonmut.genomics)

Fig4 <- rando.plot %>%
    add.cumulative.mut.layer(c.top.nonmuts,my.color="black")
ggsave("../results/gene-modules/figures/Fig4.pdf",Fig4)

## calculate more rigorous statistics than the figures.
top.nonmut.pvals <- calc.traj.pvals(gene.mutation.data,
                                    REL606.genes,
                                    unique(top.nonmut.genomics$Gene.name))
## recalculate, sampling from the same genomic regions (bins).
## results should be unchanged.
top.nonmut.pvals.loc <- calc.traj.pvals(gene.mutation.data,
                                        REL606.genes,
                                        unique(top.nonmut.genomics$Gene.name),
                                        sample.genes.by.location=TRUE)

## 2) plot top genes in hypermutators.
top.hypermut.data <- gene.mutation.data %>%
    filter(Gene %in% top.hypermut.genomics$Gene.name)

c.top.hypermut <- calc.cumulative.muts(top.hypermut.data, top.hypermut.genomics)

S1Fig <- rando.plot %>%
    add.cumulative.mut.layer(c.top.hypermut, my.color="black")
ggsave("../results/gene-modules/figures/S1Fig.pdf",S1Fig)

## calculate more rigorous statistics than the figures.
top.hypermut.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes,
                                      unique(top.hypermut.genomics$Gene.name))
## recalculate, sampling from the same genomic regions (bins).
## results should be unchanged.
top.hypermut.pvals.loc <- calc.traj.pvals(gene.mutation.data, REL606.genes,
                                          unique(top.hypermut.genomics$Gene.name),
                                          sample.genes.by.location=TRUE)

##########################################################################
## CORE GENES ANALYSIS.
##########################################################################
core.gene.assignments <- read.csv("../data/Maddamsetti2017-core-summary.csv")

core.genes <- core.gene.assignments %>% filter(panortholog == TRUE)
noncore.genes <- core.gene.assignments %>% filter(panortholog == FALSE)

core.genes.metadata.df <- filter(REL606.genes, locus_tag %in%core.genes$locus_tag)
noncore.genes.metadata.df <- filter(REL606.genes, locus_tag %in%noncore.genes$locus_tag)

core.mut.data <- filter(gene.mutation.data, locus_tag %in% core.genes$locus_tag)
noncore.mut.data <- filter(gene.mutation.data, locus_tag %in% noncore.genes$locus_tag)

core.subset.size <- length(unique(core.genes$locus_tag))
noncore.subset.size <- length(unique(noncore.genes$locus_tag))

core.rando.layer <- plot.base.layer(gene.mutation.data, REL606.genes,
                                        subset.size=core.subset.size)

c.core.muts <- calc.cumulative.muts(core.mut.data, core.genes.metadata.df)
c.noncore.muts <- calc.cumulative.muts(noncore.mut.data, core.genes.metadata.df)

coreFig <- core.rando.layer %>%
    add.cumulative.mut.layer(c.core.muts, my.color="black") %>%
    add.cumulative.mut.layer(c.noncore.muts,my.color="red")

##########################################################################
## ADHESIN GENES ANALYSIS.
##########################################################################

adhesin.genes <- REL606.genes %>%
    filter(str_detect(product, "adhesin"))

adhesin.mut.data <- filter(gene.mutation.data, locus_tag %in% adhesin.genes$locus_tag)

adhesin.subset.size <- length(unique(adhesin.genes$locus_tag))

adhesin.rando.layer <- plot.base.layer(gene.mutation.data, REL606.genes,
                                        subset.size=adhesin.subset.size)

c.adhesin.muts <- calc.cumulative.muts(adhesin.mut.data, adhesin.genes)

adhesinFig <- adhesin.rando.layer %>%
    add.cumulative.mut.layer(c.adhesin.muts, my.color="black")

##########################################################################
## GENE MODULE ANALYSIS.
##########################################################################

## First, negative results that go into the supplement.
## downstream modules in the GRN that are predictive of growth show little 
## sign of selection.

## look at accumulation of stars over time for genes in the different proteome
## sectors.
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN in each population.

## get proteome sector assignments from Hui et al. 2015 Supplementary Table 2.
## I saved a reduced version of the data.
proteome.assignments <- read.csv('../data/Hui-2015-proteome-sector-assignments.csv',
                                 as.is=TRUE) %>%
    inner_join(REL606.genes)
## add proteome assignment to gene mutation.data.
sector.mut.data <- inner_join(gene.mutation.data, proteome.assignments)
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

proteome.rando.layer <- plot.base.layer(gene.mutation.data, REL606.genes,
                                        subset.size=proteome.subset.size)

c.A.muts <- calc.cumulative.muts(A.sector.mut.data, filter(proteome.assignments,Sector.assigned=='A'))
c.S.muts <- calc.cumulative.muts(S.sector.mut.data, filter(proteome.assignments,Sector.assigned=='S'))
c.O.muts <- calc.cumulative.muts(O.sector.mut.data, filter(proteome.assignments,Sector.assigned=='O'))
c.U.muts <- calc.cumulative.muts(U.sector.mut.data,filter(proteome.assignments,Sector.assigned=='U'))
c.R.muts <- calc.cumulative.muts(R.sector.mut.data,filter(proteome.assignments,Sector.assigned=='R'))
c.C.muts <- calc.cumulative.muts(C.sector.mut.data,filter(proteome.assignments,Sector.assigned=='C'))

## plot for proteome sectors.
S2Fig <- proteome.rando.layer %>%
    add.cumulative.mut.layer(c.A.muts, my.color="black") %>%
    add.cumulative.mut.layer(c.S.muts,my.color="red") %>%
    add.cumulative.mut.layer(c.O.muts,my.color="blue") %>%
    add.cumulative.mut.layer(c.U.muts,my.color="green") %>%
    add.cumulative.mut.layer(c.R.muts,my.color="yellow") %>%
    add.cumulative.mut.layer(c.C.muts,my.color="orange")
ggsave(S2Fig, filename='../results/gene-modules/figures/S2Fig.pdf')
##########################################################################
## look at accumulation of stars over time for genes in different eigengenes
## inferred by Wytock and Motter (2018).

## get eigengene sector assignments from Wytock and Motter (2018) Supplementary File 1.
## I saved a reduced version of the data.

eigengenes <- read.csv('../data/Wytock2018-eigengenes.csv',as.is=TRUE) %>%
    inner_join(REL606.genes,eigengenes)
## add eigengene assignment to gene.mutation.data.
eigengene.mut.data <- inner_join(gene.mutation.data, eigengenes)

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

eigen.rando.layer <- plot.base.layer(gene.mutation.data,REL606.genes,
                                     subset.size=eigen.subset.size)

c.eigen1.muts <- calc.cumulative.muts(eigengene1.mut.data, filter(eigengenes,Eigengene==1))
c.eigen2.muts <- calc.cumulative.muts(eigengene2.mut.data, filter(eigengenes,Eigengene==2))
c.eigen3.muts <- calc.cumulative.muts(eigengene3.mut.data, filter(eigengenes,Eigengene==3))
c.eigen4.muts <- calc.cumulative.muts(eigengene4.mut.data, filter(eigengenes,Eigengene==4))
c.eigen5.muts <- calc.cumulative.muts(eigengene5.mut.data, filter(eigengenes,Eigengene==5))
c.eigen6.muts <- calc.cumulative.muts(eigengene6.mut.data, filter(eigengenes,Eigengene==6))
c.eigen7.muts <- calc.cumulative.muts(eigengene7.mut.data, filter(eigengenes,Eigengene==7))
c.eigen8.muts <- calc.cumulative.muts(eigengene8.mut.data, filter(eigengenes,Eigengene==8))
c.eigen9.muts <- calc.cumulative.muts(eigengene9.mut.data, filter(eigengenes,Eigengene==9))

S3Fig <- eigen.rando.layer %>%
    add.cumulative.mut.layer(c.eigen1.muts, my.color="red") %>%
    add.cumulative.mut.layer(c.eigen2.muts,my.color="orange") %>%
    add.cumulative.mut.layer(c.eigen3.muts,my.color="yellow") %>%
    add.cumulative.mut.layer(c.eigen4.muts,my.color="green") %>%
    add.cumulative.mut.layer(c.eigen5.muts,my.color="cyan") %>%
    add.cumulative.mut.layer(c.eigen6.muts,my.color="blue") %>%
    add.cumulative.mut.layer(c.eigen7.muts,my.color="violet") %>%
    add.cumulative.mut.layer(c.eigen8.muts,my.color="pink") %>%
    add.cumulative.mut.layer(c.eigen9.muts,my.color="black")
ggsave(S3Fig, filename='../results/gene-modules/figures/S3Fig.pdf')

####################################
## Now for real results.

## look at accumulation of stars over time for genes in different transcriptional
## modules inferred by Sastry et al. (2020) paper from Bernhard Palsson's group.

## Questions:
## Q1) are the regulators of the I-modulons under stronger selection than
## the genes that they regulate?

## Q2) Do any I-modulons show evidence of positive or purifying selection?
## I may have to do some kind of FDR correction for multiple hypothesis testing.

## I made this file by hand, by going through the I-modulons in imodulon_gene_names.txt,
## and the regulator annotations in modulon.pdf, both in precise-db-repo.
Imodulons.to.regulators <- read.csv("../data/rohans-I-modulons-to-regulators.csv")
Imodulon.regulators <-Imodulons.to.regulators %>%
    filter(!(is.na(regulator)))
Imodulon.regulator.mut.data <- gene.mutation.data %>%
    filter(Gene %in% Imodulon.regulators$regulator)
c.Imodulon.regulators <- calc.cumulative.muts(Imodulon.regulator.mut.data,
                                              Imodulon.regulators)

## calculate more rigorous statistics than the figures.
Imodulon.regulator.pvals <- calc.traj.pvals(gene.mutation.data, REL606.genes,
                                            unique(Imodulon.regulators$regulator))
## recalculate, sampling from the same genomic regions (bins).
## results should be unchanged.
Imodulon.regulator.pvals.loc <- calc.traj.pvals(gene.mutation.data, REL606.genes,
                                                unique(Imodulon.regulators$regulator),
                                                sample.genes.by.location=TRUE)

## Now look at genes that are regulated within Imodulons.
## I expect relaxed or purifying selection overall.

## this file was generated by reformat-I-modulons.py
genes.to.Imodulons <- read.csv("../results/gene-modules/genes-to-I-modulons.csv")
Imodulon.regulated <- genes.to.Imodulons %>% filter(!(is.na(Gene))) %>%
    ## remove any I-modulon regulators from the I-modulon genes to
    ## make an orthogonal comparison.
    filter(!(Gene %in% Imodulon.regulators$regulator))
Imodulon.regulated.mut.data <- gene.mutation.data %>%
    filter(Gene %in% Imodulon.regulated$Gene)
c.Imodulon.regulated <- calc.cumulative.muts(Imodulon.regulated.mut.data,Imodulon.regulated)

Imodulon.regulators.base.layer <- plot.base.layer(
    gene.mutation.data,
    REL606.genes,
    subset.size=length(unique(Imodulon.regulators$regulator)))

## Figure for paper:  compare I-modulon regulators to the genes they regulate.
Fig5 <- Imodulon.regulators.base.layer %>% ## null for regulators
    add.base.layer(gene.mutation.data, ## add null for regulated genes
                   REL606.genes,
                   subset.size=length(unique(Imodulon.regulated$Gene)),
                   my.color="pink") %>%
    add.cumulative.mut.layer(c.Imodulon.regulators, my.color="black") %>%
    add.cumulative.mut.layer(c.Imodulon.regulated, my.color="red")
ggsave("../results/gene-modules/figures/Fig5.pdf", Fig5)

## Make plots for each I-modulon.

## This helper function refers to several global variables.
## TODO: Best to wrap this into the context of a larger function later,
## to maintain modularity.
make.modulon.plots.helper <- function(my.I.modulon,plot.regulators=FALSE) {
    
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

    c.my.modulon.muts <- calc.cumulative.muts(my.modulon.mut.data, my.modulon.genes)
    ## for the plots, subsample based on the cardinality of the I-modulon.
    modulon.size <- length(unique(my.modulon.genes$Gene))
    modulon.base.layer <- plot.base.layer(gene.mutation.data,
                                          REL606.genes,
                                          subset.size=modulon.size,
                                          my.color="pink")
    p <- modulon.base.layer %>%
        add.cumulative.mut.layer(c.my.modulon.muts,my.color="red")

    ## if any regulators exist and flag set, then
    ## add layers for mutations in those genes.
    if ((regulator.size > 0) & plot.regulators) {
        my.regulator.mut.data <- gene.mutation.data %>%
            filter(Gene %in% my.regulators$regulator)
        c.my.regulators <- calc.cumulative.muts(my.regulator.mut.data, my.regulators)

        p <- p %>% 
            add.base.layer(gene.mutation.data, REL606.genes, my.color="grey",
                           subset.size=regulator.size) %>%
            add.cumulative.mut.layer(c.my.regulators,my.color='black')
    }
    ## add a title to the plots.
    p <- p + ggtitle(modulon.text)
    return(p)
}

## too big for PDF-- plot to separate pngs.
## to make Supplementary File S1, run "python stitchS1File.py".
png("../results/gene-modules/figures/I-modulon-plots/plot-%02d.png")
## split into groups by I-modulon, and map them to the plotting helper.
Imodulons.to.regulators %>% group_split(I.modulon) %>%
    map(.f=make.modulon.plots.helper)
dev.off()

######################################################################################
## CIS-REGULATORY EVOLUTION IN REGULATORS VERSUS EFFECTORS IN I-MODULONS.
## Let's examine cis-regulatory evolution by examining non-coding mutations.

## 1920 noncoding mutations with a gene annotation.
noncoding.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='noncoding') %>%
    filter(Gene != 'intergenic')

## Let's split into I-modulon regulators and regulated genes.
## 48 mutations associated with 70 I-modulon regulators.
Imodulon.regulator.noncoding.data <- noncoding.mutation.data %>%
    filter(Gene %in% Imodulon.regulators$regulator)
## 542 associated with 1394 I-modulon regulated genes.
Imodulon.regulated.noncoding.data <- noncoding.mutation.data %>%
    filter(Gene %in% Imodulon.regulated$Gene)
## so seems like I-modulon regulators have more noncoding hits than expected.
binom.test(x=48, n=(48+542),p=(70/(70+1394)))

## 415 with greater than 1 hit.
noncoding.by.gene <- noncoding.mutation.data %>%
    group_by(Gene) %>% summarize(count=n()) %>% arrange(desc(count)) %>%
    left_join(REL606.genes) %>%
    filter(count>1)

## 134 with greater than 2 hits.
noncoding.by.gene2 <- noncoding.by.gene %>% filter(count>2)

## 44 with greater than 3 hits.
noncoding.by.gene3 <- noncoding.by.gene %>% filter(count>3)

## 10 I-modulon regulators with multiple non-coding hits
Imodulon.regulator.noncoding.data2 <- noncoding.by.gene %>%
    filter(Gene %in% Imodulon.regulators$regulator)

## 115 I-modulon.regulated with multiple non-coding hits.
Imodulon.regulated.noncoding.data2 <- noncoding.by.gene %>%
    filter(Gene %in% Imodulon.regulated$Gene)

## what are the multiple non-coding hits that are not in an I-modulon?
## 290 of these!
non.Imodulon.noncoding.hits <- noncoding.by.gene %>%
    filter(!(Gene %in% Imodulon.regulators$regulator)) %>%
    filter(!(Gene %in% Imodulon.regulated$Gene))

## 94 of these.
non.Imodulon.noncoding.hits2 <- non.Imodulon.noncoding.hits %>%
    filter(count>2)             
##########################################################################
## PSEUDOMONAS HYPERMUTATOR EXPERIMENTAL EVOLUTION

## As an additional test, let's see if regulators are under selection in the
## hypermutator populations of Pseudomonas aeruginosa
## evolving in the lab under antibiotic selection.

## Data comes from Mehta et al. (2018):
## The Essential Role of Hypermutation in Rapid Adaptation to Antibiotic Stress.

PAO11.genes <- read.csv("../results/PAO11_IDS.csv")
## let's find all genes that match these keywords:
## regulator, regulatory, regulation.
PAO11.regulators <- PAO11.genes %>%
    filter(str_detect(product, 'regulator|regulatory|regulation'))
## OK. there are 424 such genes.

## Now, let's run STIMS on the Mehta data. Ignore control pop.
Mehta.data <- read.csv("../results/Mehta2018-hypermutators.csv") %>%
    ## remove intergenic mutations.
    filter(!str_detect(Gene,'/')) %>%
    left_join(PAO11.genes) %>%
    mutate(Day = t0) %>%
    filter(Population != 'control') %>%
    mutate(Population = factor(Population)) %>%
    mutate(Population = recode(Population, `replicate1` = "Replicate 1",
                               `replicate2` = "Replicate 2")) %>%
    ## filter any rows with NAs.
    drop_na()

## show both replicate populations in Supplementary Figure S4.
Mehta.regulator.mut.data <- Mehta.data %>%
    filter(Gene %in% PAO11.regulators$Gene)
## 188 mutations occurring in annotated regulators.

c.Mehta.regulators <- calc.Mehta.cumulative.muts(Mehta.regulator.mut.data, PAO11.regulators)
Mehta.base.layer <- plot.Mehta.base.layer(Mehta.data, PAO11.genes,
                                          subset.size=nrow(PAO11.regulators))
S4Fig <- Mehta.base.layer %>%
    add.Mehta.cumulative.mut.layer(c.Mehta.regulators, my.color="black")
ggsave("../results/gene-modules/figures/S4Fig.pdf", S4Fig, width=5, height=5)

## combine both replicates for Figure 6.
combined.Mehta.data <- Mehta.data %>%
    mutate(Population = recode(Population,
                               `Replicate 1` = "Replicate populations combined",
                               `Replicate 2` = "Replicate populations combined"))

combined.Mehta.regulator.mut.data <- combined.Mehta.data %>%
    filter(Gene %in% PAO11.regulators$Gene)

c.combined.Mehta.regulators <- calc.Mehta.cumulative.muts(
    combined.Mehta.regulator.mut.data, PAO11.regulators)

combined.Mehta.base.layer <- plot.Mehta.base.layer(combined.Mehta.data,
                                                   PAO11.genes,
                                                   subset.size=nrow(PAO11.regulators))
Fig6 <- combined.Mehta.base.layer %>%
    add.Mehta.cumulative.mut.layer(c.combined.Mehta.regulators, my.color="black")
ggsave("../results/gene-modules/figures/Fig6.pdf", Fig6, width=5, height=5)

## calculate rigorous statistics for Figure 6.
Mehta.pvals <- calc.Mehta.traj.pvals(combined.Mehta.data, PAO11.genes,
                                     unique(PAO11.regulators$Gene))

##########################################################################
## Examples plots for the new Figure 1.
## The new figure 1 shows:
## 1) the nature of the input data,
## 2) the test statistic
## 3) what comparisons are being performed.
## It does so by showing:
## A) the STIMS pipeline,
## B) the difference between black lines and gray area,
## C) the probability under the null statistical model.

## 71 ribosome-associated proteins. I found these by searching the REL606 genbank file
## using the keyword "ribosom" to match "ribosomal" and "ribosome".
## There's 72 by manual search, but 71 in REL606.genes.

ribosome.associated.prot.genes <- REL606.genes %>%
    filter(str_detect(product, "ribosom"))

ribosome.associated.prot.subset.size <- length(ribosome.associated.prot.genes$locus_tag)

## for simplicity, just plot data for Ara-1.
m1.pop.gene.mutation.data <- filter(gene.mutation.data, Population == "Ara-1") %>%
    ## change levels of the population factor so that only Ara-1 is represented.
    mutate(Population = factor(
               Population,
               levels = unique(Population)))

m1.ribosome.associated.prot.mut.data <- filter(
    m1.pop.gene.mutation.data,
    locus_tag %in% ribosome.associated.prot.genes$locus_tag)

c.m1.ribosome.associated.prot.muts <- calc.cumulative.muts(
    m1.ribosome.associated.prot.mut.data,
    ribosome.associated.prot.genes,
    reset.pop.levels = FALSE)

## to make the plots uniform, plot the rando layer as the invisible background.
m1.blank.rando.layer <- plot.base.layer(m1.pop.gene.mutation.data, REL606.genes,
                                        subset.size=ribosome.associated.prot.subset.size,
                                        my.color="white")

m1.ribosome.associated.prot.Fig <- m1.blank.rando.layer %>%
    add.cumulative.mut.layer(c.m1.ribosome.associated.prot.muts, my.color="black")
ggsave("../results/gene-modules/figures/Fig0-panels/Fig0-STIMS-1.pdf",
       m1.ribosome.associated.prot.Fig,height=4.2, width=4.2)

## plot this panel to add to the figure.
m1.rando.layer <- plot.base.layer(m1.pop.gene.mutation.data, REL606.genes,
                                  subset.size=ribosome.associated.prot.subset.size)
ggsave("../results/gene-modules/figures/Fig0-panels/Fig0-STIMS-2.pdf",
       m1.rando.layer,height=4.2, width=4.2)


exampleFig <- m1.rando.layer %>%
    add.cumulative.mut.layer(c.m1.ribosome.associated.prot.muts, my.color="black")
ggsave("../results/gene-modules/figures/Fig0-panels/Fig0-STIMS-3.pdf",
       exampleFig,height=4.2, width=4.2)

