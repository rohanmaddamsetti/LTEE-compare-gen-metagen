## generate-SLiM-genome-files.R by Nkrumah Grant and Rohan Maddamsetti
## Create a simple genome.csv file akin to REL606_IDs.csv,
## and create files for each module.

## Each of the genomic elements is divided into 1000 genes, each of size 100.
n <- 400
prefix <- "g"
suffix <- seq(1:n)
Gene <- paste(prefix, suffix, sep="")

end <- c()
start <- c()

for (i in 1:n) {
    end <- c(end, (i*1000))
    start <- c(end - 1000)
    gene_length <- c(end - start)
}

## There are 4 genomic elements in the SLiM genome. Each of them are 100000 basepairs in
## length and are adjacent to each other with no spaces in between, so gene 1
## would fall between position 0 and 100000, and gene 2 between 100001 and
## 200001, and so forth. This is the structure we spoke about previously. My
## thought on this is that STIMS would sample a random length of positions in
## each genomic element we are testing (e.g. one under positive selection) and
## compare that to a null distribution created by bootstrapping across genomic
## elements sampled of equal length. I have added a “Gene” column to the dataframe
## for each of the genomic elements. See below the key for what each genomic
## element represents, and the relative ratio of mutation effects, which is also
## included in my SLiM Eidos script. 

## neutral, beneficial, deleterious 
## g1 = 1, 0, 0, 
## g2 = 0.85, 0.10, 0.05 
## g3 = 0.45, 0.40, 0.05 
## g4 = 0.50, 0,  0.50

SLiM.genes <- data.frame(Gene, start, end, gene_length) %>%
    ## Let's annotate the genome file.
    mutate(Module = c(rep("Neutral", 100), rep("Weak_Positive",100),
                      rep("Positive", 100), rep("Purifying",100))) %>%
    mutate(PercentageNeutral = c(rep(1, 100), rep(0.85, 100),
                                 rep(0.45, 100), rep(0.50,100))) %>%
    mutate(PercentageBeneficial = c(rep(0, 100), rep(0.1, 100),
                                    rep(0.4, 100), rep(0,100))) %>%
    mutate(PercentageDeleterious = c(rep(0, 100), rep(0.05, 100),
                                     rep(0.05, 100), rep(0.50,100)))

write.csv(SLiM.genes, "../results/SLiM-results/SLiM_geneIDs.csv",
          quote = F, row.names = F)

neutral.module <- filter(SLiM.genes, Module == "Neutral")
weak.positive.selection.module <- filter(SLiM.genes, Module == "Weak_Positive")
positive.selection.module <- filter(SLiM.genes, Module == "Positive")
purifying.selection.module <- filter(SLiM.genes, Module == "Purifying")
## write the modules to file.
write.csv(file="../results/SLiM-results/SLiM_neutral_module.csv",
          neutral.module)
write.csv(file="../results/SLiM-results/SLiM_weak_positive_module.csv",
          weak.positive.selection.module)
write.csv(file="../results/SLiM-results/SLiM_positive_module.csv",
          positive.selection.module)
write.csv(file="../results/SLiM-results/SLiM_purifying_module.csv",
          purifying.selection.module)
