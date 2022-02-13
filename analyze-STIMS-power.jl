"""
analyze-STIMS-power.jl by Rohan Maddamsetti.

What variables affect STIMS statistical power to detect selection? 
To answer this question, we asked how the following variables affect STIMS Type I error 
(detecting selection when no such pattern exists) and Type II error rate 
(ability to detect selection given the pattern exists).

1) length of the time series. Vary to mimic the kind of experiment people could do: 
100, 500, 1000, 2000, 5000 generations. No need to run multiple simulations, 
just cut the data at different timepoints.
2) sampling interval. 10 gen, 50 gen, 100 gen, 500 gen. Again, no need to run multiple 
simulations, since we can cut the data at different timepoints.
3) number of genes chosen from the ground truth module: 
1 gene, 10 genes, 25 genes, 100 genes.
4) metagenomic filtering threshold: no filtering, 0.1% 1%, 10%. 

Fixed parameters for simulations:
Fix populations at Ne = 10^6
Fix 2 different mutation rates: mu = 10^-8 and 10^-10.
Fix the genome design, and always run STIMS on beneficial module, 
neutral module, and deleterious module.

Fixed parameters when the other parameters are varied:
Sampling interval = 100 generations.
Time Series = 5,000 generations
Number of genes chosen from ground truth module = 100.
Metagenomic filtering threshold = 1%

3 modules X 2 mutation rates X 4 experimental variables.
So, I will make 4 sets of 6 figure panels.

Program design:

For each replicate dataset:
-- Take SLiM data.
-- Downsample by sampling interval.
-- Filter by metagenomic sampling threshold.
-- Filter the total length of the time series.
-- Put into the STIMS data format.
-- For each module:
      -- Take the first n genes of the module.
      -- Run STIMS.jl.
      -- Get the p-value, and add one row to the dataframe.
-- Use the dataframe to make figures using ggplot2.

"""
