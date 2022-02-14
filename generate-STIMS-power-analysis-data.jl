"""
generate-STIMS-power-analysis-data.jl by Rohan Maddamsetti.

What variables affect STIMS statistical power to detect selection? 
To answer this question, we asked how the following variables affect STIMS Type I error 
(detecting selection when no such pattern exists) and Type II error rate 
(ability to detect selection given the pattern exists).

1) length of the time series. Vary to mimic the kind of experiment people could do: 
500, 1000, 2000, 5000 generations. No need to run multiple simulations, 
just cut the data at different timepoints.
2) sampling interval. 25 gen, 50 gen, 100 gen, 500 gen. Again, no need to run multiple 
simulations, since we can cut the data at different timepoints.
3) number of genes chosen from the ground truth module: 
1 gene, 10 genes, 25 genes, 50, 100 genes.
4) metagenomic filtering threshold: no filtering, 0.1%, 1%, 5%, 10%. 

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
-- write the dataframe to file.

Then, use the R script analyze-STIMS-power.R to make figures using ggplot2.
"""

include("STIMS.jl")
using DataFrames, DataFramesMeta, CSV, RCall


## IMPORTANT: we are assuming that each gene is 1000 bp long.
function annotate_Gene_per_mutation(df; gene_length = 1000)
    ## annotate the Gene for each mutation.
    annotated_df = @chain df begin
        @rtransform(:GeneBin = div(:Position, gene_length) + 1)
        @rtransform(:GeneBin = "g"*string(:GeneBin)) 
        rename(:GeneBin => :Gene)
    end
end


## This function nicely formats all mutations over time into a dataframe.
## Further processing is needed to generate STIMS input.
function SLiM_output_to_dataframe(SLiM_file, pop_name, timeseries_length,
                                  sampling_interval, freq_threshold, Ne)

    SLiM_output = CSV.read(SLiM_file, DataFrame; header = false, delim=" ")

    SLiM_df = @chain SLiM_output begin
        ## Remove unnecessary columns from SLiM output. 
        @select(:Column2, :Column5, :Column6, :Column7, :Column10, :Column11, :Column12)
        rename(:Column2 => :Generation)
        rename(:Column5 => :ID)
        rename(:Column6 => :Annotation)
        rename(:Column7 => :Position)
        rename(:Column10 => :Population)
        rename(:Column11 => :t0)
        rename(:Column12 => :prevalence)
        @transform(:Annotation = replace(:Annotation,
                                         "m1" => "beneficial",
                                         "m2" => "deleterious",
                                         "m3" => "neutral",
                                         "m4" => "background"))
        @transform(:Population = replace(:Population, "p1" => pop_name))
        ## annotate the Gene for each mutation.
        annotate_Gene_per_mutation()
        ## Convert prevalence to allele frequency.
        @rtransform(:allele_freq = :prevalence/Ne)
        ## Filter mutations with allele frequencies above the sampling threshold.
        @rsubset(:allele_freq >= freq_threshold)
        ## Filter mutations within the given timeseries length.
        @rsubset(:Generation <= timeseries_length)
        ## Filter mutations per the sampling interval.
        @rsubset(mod(:Generation, sampling_interval) == 0)
        sort([:ID, :Generation]) ## arrange each mutation by generation
    end
end


## This function makes STIMS input-- only one mutation per row,
## representing the first time the mutation was observed in the population.
## set the default parameters that are held fixed when one parameter is varied.
function SLiM_output_to_STIMS_input(SLiM_file, pop_name;
                                    timeseries_length = 5000,
                                    sampling_interval = 100,
                                    freq_threshold = 0.01,
                                    Ne = 1e6)
    
    SLiM_df = SLiM_output_to_dataframe(SLiM_file, pop_name, timeseries_length,
                                       sampling_interval, freq_threshold, Ne)

    STIMS_input_df = @chain SLiM_df begin
        ## IMPORTANT: There can only be one row per mutation.
        groupby([:ID, :Annotation, :Position, :Gene, :Population, :t0])
        combine(nrow)
        @select(:ID, :Annotation, :Position, :Gene, :Population, :t0)
        sort(:ID) ## to double-check that each ID is unique
    end
    return(STIMS_input_df)
end


function getMutationRate(SLiM_file)
    fname = basename(SLiM_file)
    mutation_rate_field = split(fname, "_")[3]
    logMu = split(mutation_rate_field, "-")[2]
    if (logMu == "10")
        pop_name = "Nonmutator"
    elseif (logMu == "8")
        pop_name = "Hypermutator"
    else
        pop_name = "Unknown_mu_rate"
    end
    return pop_name
end


function GeneratePowerAnalysisDataFrame()
    ## Add rows to this dataframe, one by one.
    power_analysis_df = DataFrame(MutatorStatus = String[],
                                  GeneModule = String[],
                                  Variable = String[],
                                  VariableValue = Float64[],
                                  Count = Float64[],
                                  Pval = Float64[])
    ## Variable parameters to test.
    timeseries_length_vec = [500, 1000, 2000, 5000]
    sampling_interval_vec = [25, 50, 100, 500]
    num_genes_from_module_vec = [1, 10, 25, 50, 100]
    metagenomic_filtering_threshold_vec = [0, 0.001, 0.01, 0.05, 0.1]
    
    outdir = "../results/SLiM-results"
    genome_metadata_csv_path = joinpath(outdir, "SLiM_geneIDs.csv")
    genome_metadata = CSV.read(genome_metadata_csv_path, DataFrame)
    
    gene_set_path_dict = Dict("positive" => joinpath(outdir, "SLiM_positive_module.csv"),
                              "neutral" => joinpath(outdir, "SLiM_neutral_module.csv"),
                              "purifying" => joinpath(outdir, "SLiM_purifying_module.csv"))
    gene_module_df_dict = Dict(k => CSV.read(v, DataFrame) for (k,v) in gene_set_path_dict)
    
    SLiM_dataset_dir = "../results/SLiM-results/SLiM-replicate-runs"
    SLiM_dataset_vec = readdir(SLiM_dataset_dir; join=true)
    
    for SLiM_file in SLiM_dataset_vec
        ## name as "Hypermutator" or "Nonmutator" based on the filename.
        pop_name = getMutationRate(SLiM_file)
        
        for exp_length in timeseries_length_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; timeseries_length = exp_length)
            for (module_type, module_df) in gene_module_df_dict         
                pval_df = STIMS.RunSTIMS_on_data(STIMS_input, genome_metadata, module_df)
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, "timeseries_length",
                                          exp_length, pval_df.count[1], pval_df.pval[1]])
            end
        end
        
        for sampling_int in sampling_interval_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; sampling_interval = sampling_int)
            for (module_type, module_df) in gene_module_df_dict
                pval_df = STIMS.RunSTIMS_on_data(STIMS_input, genome_metadata, module_df)
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, "sampling_interval",
                                          sampling_int, pval_df.count[1], pval_df.pval[1]])
            end
        end
        
        for threshold in metagenomic_filtering_threshold_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; freq_threshold = threshold)
            for (module_type, module_df) in gene_module_df_dict
                pval_df = STIMS.RunSTIMS_on_data(STIMS_input, genome_metadata, module_df)
                println(pval_df)
                push!(power_analysis_df, [pop_name,  module_type, "filtering_threshold",
                                          threshold, pval_df.count[1], pval_df.pval[1]])
            end
        end
        ## Note that this case is a bit different.
        for num_genes in num_genes_from_module_vec
            STIMS_input = SLiM_output_to_STIMS_input(SLiM_file, pop_name)
            for (module_type, module_df) in gene_module_df_dict
                ## Only run on the first num_genes in the module.
                filtered_gene_module_df = module_df[1:num_genes, :]
                pval_df = STIMS.RunSTIMS_on_data(STIMS_input, genome_metadata,
                                                 filtered_gene_module_df)
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, "num_genes_from_module",
                                          num_genes, pval_df.count[1], pval_df.pval[1]])
            end
        end
    end
    return(power_analysis_df)
end


## make the DataFrame for the power analysis
power_analysis_df = GeneratePowerAnalysisDataFrame()
## and write to file. Make figures with analyze-STIMS-power.R.
CSV.write("../results/SLiM-results/STIMS-power-analysis-data.csv", power_analysis_df)
