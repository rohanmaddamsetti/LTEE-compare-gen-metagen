## SLiM-plot-mutation-trajectories.R
## Nkrumah Grant and Rohan Maddamsetti
## Take SLiM output data and create a mutation trajectory plot.


library(tidyverse)
library(gghighlight)
library(cowplot)
library(magick)
library(ggplotify)

######Chunk 1######
#Set working directory to the location of SLiM output.
setwd("/Users/nkrumahgrant/Desktop")


#Read files into the global environment.
#For files with fixed mutations, remove the first two rows from the dataframe, which just
#corresponds to the file name and the generation of SLiM output. All SliMulations were 
#completed over of 10,000 generations. 

d1.1 <- read.csv("20211211.Nonmutator.NeutralT.txt", header = F, sep = " ")
d1.2 <- read.csv("20211211.Nonmutator-NuetralTFixed.txt", header = F, sep = " ", skip = 2)

d2.1 <- read.csv("20211211.Hypermutator.NeutralT.txt", header = F, sep = " ") 
d2.2 <- read.csv("20211211.Hypermutator-NuetralTFixed.txt", header = F, sep = " ", skip = 2)

#Remove unnecessary columns from SLiM output. 
#Nonmutator
d1.1 <- d1.1 %>% select (c("V2", "V5" , "V6", "V7", "V10", "V11", "V12"))
d1.2 <- d1.2 %>% select (c("V2", "V4", "V7", "V8", "V9"))

#Hypermutator
d2.1 <- d2.1 %>% select (c("V2", "V5" , "V6", "V7", "V10", "V11", "V12"))
d2.2 <- d2.2 %>% select (c("V2", "V4", "V7", "V8", "V9"))

#Rename columns for all dataframes.
#X.1 reflects all mutations that were segregating in the population as were sampled every 100 
#generations. X.2 reflects all mutations that fixed in the population. 
#Where X = 1|Nonmutator and X = 2|Hypermutator scenarios. 
names(d1.1) <-  (c("Generation", "ID", "Annotation", "Position", "Population", "t0", "prevalance"))
names(d2.1) <-  (c("Generation", "ID", "Annotation", "Position", "Population", "t0", "prevalance"))

names(d1.2) <- (c("ID", "Position", "Population", "t0", "MutFixed"))
names(d2.2) <- (c("ID", "Position", "Population", "t0", "MutFixed"))

#Convert position column to factor
#Nonmutator
d1.1$Position <- as.factor(d1.1$Position)
d1.2$Position <- as.factor(d1.2$Position)

#Hypermutator
d2.1$Position <- as.factor(d2.1$Position)
d2.2$Position <- as.factor(d2.2$Position)

#rekey factor level names of mutation identifiers m1, m2, m3 to neutral, beneficial, 
#and deleterious, respectively. 
#Nonmutator 
levels(d1.1$Annotation) <- list(neutral = "m1", beneficial = "m2", deleterious = "m3")
levels(d1.1$Population) <- "Nonmutator"

#Hypermutator
levels(d2.1$Annotation) <- list(neutral = "m1", beneficial = "m2", deleterious = "m3")
levels(d2.1$Population) <- "Hypermutator"

#Join datasets

#Nonmutator 
d3.1 <- full_join(d1.1, d1.2)

#Hyperutator 
d4.1 <- full_join(d2.1, d2.2)

#Joining the two dataframes created columns with NA values. This is because there are some 
#uncommon columns between the tew dataframes. Here, I manually correct that. 

#Nonmutator
d3.1$Annotation[is.na(d3.1$Annotation)]<- as.factor("beneficial")
d3.1$Generation[is.na(d3.1$Generation)]<- as.integer("10000")
d3.1$prevalance[is.na(d3.1$prevalance)]<- as.integer("500000")
d3.1$MutFixed[is.na(d3.1$MutFixed)] <- as.integer("0")

#Hypermutator
d4.1$Annotation[is.na(d4.1$Annotation)]<- as.factor("beneficial")
d4.1$Generation[is.na(d4.1$Generation)]<- as.integer("10000")
d4.1$prevalance[is.na(d4.1$prevalance)]<- as.integer("500000")
d4.1$MutFixed[is.na(d4.1$MutFixed)] <- as.integer("0")

#Create new column that converts prevalance to allele frequency
#This will be used to generate mutation trajectories. These figures will be placed in 
#supplementary material. The denominator is the number of individuals in the SliMulation. 

#Nonmutator 
d3.1 <- d3.1 %>% mutate("all_freq" = prevalance/5e5) 
#Hypermutator
d4.1 <- d4.1 %>% mutate("all_freq" = prevalance/5e5) 

#Filter mutations with allele frequencies above 5% threshold. 
#This may be another parameter to look at. If methods are developed
#where every single mutation that ever occurs in a population
#can be identified, then STIMS will fail to detect selection at high
#mutation rates. This may even occur at ancestral mutation rates, too. 
#Maybe look at 5% and 10% frequency thresholds. 

#Nonmutator
d3.1 <- d3.1 %>% filter(all_freq > 0.05)

#Hypermutator
d4.1 <- d4.1 %>% filter(all_freq > 0.05)

#Generate STIMs input dataframe 
#tf is a dummy variable. Sampling time every 100 generations. A mutation can fix or 
#go extint at anytime in that interval.

#Nonmutator
d3.2 <- as.data.frame(d3.1 %>% group_by(Population,ID,Annotation,Position,t0) %>% 
  summarise(n = n()) %>% 
  mutate(tf = t0 + 100 * n) %>% 
  select(-c(n)))

d3.2$Position <- as.character(d3.2$Position)

d3.3 <- d3.2 %>%
  filter(Population != "p1") %>% 
  mutate(P1000 = trunc(as.numeric(Position)/1000)+1) %>% 
  mutate(Gene = paste("g", P1000, sep = "")) %>% 
  select(-c(P1000))

#This is the nonmutator dataset that will be used as input for STIMS. 
#This file should be placed in the "results" directory. 
#SubNeutral in the file name relects the fact that synonymous mutations are being converted to 
#substitution objects in the simulation. Doing so simulates hitchiking dynamics. 
write.csv(d3.3, "SLIM-10000gen-FivePercent-Nonmut-SubNeutral.csv", quote= F, row.names = F)

#Hypermutator 
d4.2 <- as.data.frame(d4.1 %>% group_by(Population,ID,Annotation,Position,t0) %>% 
                        summarise(n = n()) %>% 
                        mutate(tf = t0 + 100 * n) %>% 
                        select(-c(n)))

d4.2$Position <- as.character(d4.2$Position)

d4.3 <- d4.2 %>%
  filter(Population != "p1") %>% 
  mutate(P1000 = trunc(as.numeric(Position)/1000)+1) %>% 
  mutate(Gene = paste("g", P1000, sep = "")) %>% 
  select(-c(P1000))

#This is the hypermutator dataset that will be used as input for STIMS. 
#This file should be placed in the "results" directory. 
write.csv(d4.3, "SLIM-10000gen-FivePercent-Hypermut-SubNeutral.csv", quote= F, row.names = F)


######Chunk 2######
#Make mutation trajectory plots

#Nonmutator
p1 <- d3.1 %>% 
  filter(Population != "p1") %>% 
  ggplot(aes(Generation, all_freq, colour = as.factor(ID)))+
  geom_line() +
  gghighlight(max(all_freq) > .90, max(as.numeric(Annotation)) > 0) +  
  facet_wrap(~Population, 2,1, labeller = label_parsed) +
  xlab("Generation, t") +
  ylab("Allele frequency f(t)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 14))

#Hypermutator
p2 <- d4.1 %>% 
  filter(Population != "p1") %>% 
  ggplot(aes(Generation, all_freq, colour = as.factor(ID)))+
  geom_line() +
  gghighlight(max(all_freq) > .90, max(as.integer(Annotation)) == 2) + 
  #facet_wrap(~ID + Population, 2,1) +
  facet_wrap(~Population, 2,1, labeller = label_parsed) +
  xlab("Generation, t") +
  ylab("Allele frequency f(t)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_blank()) +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 14))


#https://wilkelab.org/cowplot/articles/plot_grid.html
p1 <- plot_grid(p1, labels = 'A', label_size = 14, ncol=1)
p2 <- plot_grid(p2, labels = 'B', label_size = 14, ncol=1)

NMneu <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Nonmutator-neutral.jpeg"))
NMwp <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Nonmutator-weak-positive.jpeg"))
NMpos <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Nonmutator-positive.jpeg"))
NMpur <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Nonmutator-purifying.jpeg"))
p4 <- plot_grid(NMneu,NMwp, NMpos, NMpur, labels = c('C', 'D', 'E', 'F'), label_size = 14, ncol=2)

HMneu <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Hypermutator-neutral.jpeg"))
HMwp <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Hypermutator-weak-positive.jpeg"))
HMpos <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Hypermutator-positive.jpeg"))
HMpur <- as.ggplot(image_read("/Users/nkrumahgrant/Desktop/Colloborations/RohanCollaborations/PaperTwo/STIMS-Simulations/results/gene-modules/STIMS-jl-test-figures/STIMS-Hypermutator-purifying.jpeg"))
p5 <- plot_grid(HMneu,HMwp, HMpos, HMpur, labels = c('G', 'H', 'I', 'J'), label_size = 14, ncol=2)

plot_grid(p1, p4, p2, p5, ncol = 2)

