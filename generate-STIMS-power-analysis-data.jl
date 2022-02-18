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
      -- Run STIMS.
      -- Get the p-value, and add one row to the dataframe.
-- write the dataframe to file.

On HPC, use the following command to run:
sbatch -t 96:00:00 -c 16 --mem-per-cpu=32G --wrap="julia --threads 16 generate-STIMS-power-analysis-data.jl 16"

Then, use the R script analyze-STIMS-power.R to make figures using ggplot2.
"""

using DataFrames, DataFramesMeta, CSV, StatsBase, FLoops, ArgParse

################################################################################
## FUNCTIONS FROM STIMS.jl. I copied these functions into this script to remove
## the problematic RCall dependency when running on the Duke Compute Cluster.

function make_cumulative_muts_pop_df(gene_module_mutations_df, pop, final_time)
    gene_module_mutations_in_pop_df = @rsubset(gene_module_mutations_df,
                                               :Population == pop)
    
    if (nrow(gene_module_mutations_in_pop_df) == 0)
        ## no mutations in this pop
        almost_done_df = DataFrame(Population = pop,
                                   t0 = final_time,
                                   count = 0,
                                   cs = 0)
    else
        summary_df = @chain gene_module_mutations_in_pop_df begin
            groupby([:Population, :t0])
            combine(nrow => :count)
            sort(:t0)
            transform(:count => cumsum => :cs)
        end
        
        ## since the final_time is not in ret_df, add one final row,
        ## for nicer plots.
        final_row_df = DataFrame(Population = pop,
                                 t0 = final_time,
                                 count = 0,
                                 cs = maximum(summary_df.cs))
        
        almost_done_df = vcat(summary_df, final_row_df)
    end
    ## add an row for Time == 0 (for nicer plots).
    init_row_df = DataFrame(Population = pop,
                            t0 = 0,
                            count = 0,
                            cs = 0)
    
    pop_df = vcat(init_row_df, almost_done_df)
    return pop_df
end


function calc_cumulative_muts(gene_module_df,
                              gene_mutation_data, genome_metadata, pop_level_vec)
    ## look at accumulation of stars over time
    ## in other words, look at the rates at which the mutations occur over time.
    ## To normalize, we need to supply the number of sites at risk
    ## (such as sum of gene length).

    gene_module_mut_data = @rsubset(gene_mutation_data, :Gene in gene_module_df.Gene)
    gene_module_metadata = @rsubset(genome_metadata, :Gene in gene_module_df.Gene)
    
    ## normalize by the total length of genes
    ## in the given module (in d.metadata).
    my_genes = @chain gene_module_metadata begin
        select(:Gene, :gene_length)
        unique([:Gene, :gene_length])
    end
    
    normalization_constant = sum(my_genes.gene_length)
    
    final_time = maximum(gene_mutation_data.t0) + 1
    ## + 1 so that final_time is outside of the actual data.
    
    ## we are going to concatenate the summarized data to this empty DataFrame.
    c_dat = DataFrame()
    
    for pop in pop_level_vec
        pop_df = make_cumulative_muts_pop_df(gene_module_mut_data, pop, final_time)
        ## concatenate the pop_df to cumulative_df.
        c_dat = vcat(c_dat, pop_df)
    end
    
    transform!(c_dat, :cs => ByRow(cs -> cs/normalization_constant) => :normalized_cs)
    ## possible TODO: remove any NA values that arise?
    
    return c_dat
end


function generate_cumulative_mut_subset(gene_mutation_data, genome_metadata,
                                        pop_level_vec,
                                        subset_size)
    ## This function calculates cumulative mutations for a random gene set.
    rando_genes = sample(genome_metadata.Gene, subset_size; replace=false)
    rando_genes_df = @rsubset(genome_metadata, :Gene in rando_genes)
    
    c_mut_subset = calc_cumulative_muts(rando_genes_df, gene_mutation_data, genome_metadata, pop_level_vec)
    return c_mut_subset
end


function calc_traj_pvals(gene_module_df,
                         gene_mutation_data, genome_metadata,
                         pop_level_vec; N = 10000, ncores = 4)    
    #= calculate the tail probabilities of the true cumulative mutation trajectory
    of a given vector of genes (a 'module'), based on resampling
    random sets of genes. Returns the upper tail of null distribution,
    or P(random trajectory >= the actual trajectory).
    Output: a dataframe with three columns: Population, count, p.val.
    =#

    ## each sample has the same cardinality as gene_module_df.Gene.
    subset_size = length(gene_module_df.Gene)

    data_trajectory = calc_cumulative_muts(gene_module_df,
                                           gene_mutation_data,
                                           genome_metadata,
                                           pop_level_vec)
    
    data_final_normalized_cs_summary = @chain data_trajectory begin
        groupby(:Population)
        combine(:normalized_cs => maximum => :data_final_normalized_cs)
    end

    @floop ThreadedEx(basesize = N ÷ ncores) for _ in 1:N 
        randomized_trajectory = generate_cumulative_mut_subset(gene_mutation_data,
                                                               genome_metadata,
                                                               pop_level_vec,
                                                               subset_size)
        
        randomized_final_normalized_cs_summary = @chain randomized_trajectory begin
            groupby(:Population)
            combine(:normalized_cs => maximum => :randomized_final_normalized_cs)
            ## compare the randomized subset to the actual data:
            innerjoin(data_final_normalized_cs_summary, on = :Population)
            combine(:, [:randomized_final_normalized_cs, :data_final_normalized_cs] =>
                    ByRow((a, b) -> a > b) => :greater_than_data)
        end
        
        greater_than_data_vec = randomized_final_normalized_cs_summary.greater_than_data
        @reduce() do (count_vec = zeros(length(pop_level_vec)); greater_than_data_vec)
            ## now update the counts.
            count_vec .+= greater_than_data_vec
        end
    end
    
    uppertail_prob_df = DataFrame("Population" => pop_level_vec,
                                  "count" => count_vec,
                                  "pval" => count_vec/N); 
    return uppertail_prob_df
end


function get_middle_trajectories(bootstrapped_trajectories, alphaval, N)
    ## filter out the top alphaval/2 and bottom alphaval/2 trajectories
    ## from each population, for a two-sided test.
    ## usual default is alphaval == 0.05.
    
    trajectory_summary = @chain bootstrapped_trajectories begin
        groupby([:bootstrap_replicate, :Population])
        combine(:normalized_cs => maximum => :final_norm_cs)
    end

    num_traj_to_slice = Integer(N * alphaval) ÷ 2
    
    top_trajectories = @chain trajectory_summary begin
        sort([:Population, :final_norm_cs])
        groupby(:Population)
        combine(x -> last(x, num_traj_to_slice))
        select(:Population, :bootstrap_replicate)
        @rtransform(:in_top = true)
    end

    bottom_trajectories = @chain trajectory_summary begin
        sort([:Population, :final_norm_cs])
        groupby(:Population)
        combine(x -> first(x, num_traj_to_slice))
        select(:Population, :bootstrap_replicate)
        @rtransform(:in_bottom = true)
    end

    middle_trajectories = @chain bootstrapped_trajectories begin
        leftjoin(top_trajectories, on = [:Population, :bootstrap_replicate])
        leftjoin(bottom_trajectories, on = [:Population, :bootstrap_replicate])
        @rsubset(ismissing(:in_top))
        @rsubset(ismissing(:in_bottom))
        select(Not(:in_top))
        select(Not(:in_bottom))
    end
    
    return middle_trajectories
end


## take dataframes directly as input, and don't make the plot.
function RunSTIMS_on_data(mutation_data, genome_metadata, gene_module_df, ncores = 4)
   
    gene_mutation_data = @chain mutation_data begin
        transform(:t0 => ByRow(t -> t/1000) => :Time)
        innerjoin(genome_metadata, on = :Gene)
        @rsubset(:Gene != "intergenic")
    end
    
    pop_level_vec = unique(gene_mutation_data.Population)
    
    ## make sure that only loci in genome_metadata are analyzed.
    @rsubset!(gene_module_df, :Gene in genome_metadata.Gene)

    c_gene_module = calc_cumulative_muts(gene_module_df,
                                         gene_mutation_data,
                                         genome_metadata,
                                         pop_level_vec)

    pvals = calc_traj_pvals(gene_module_df,
                            gene_mutation_data,
                            genome_metadata,
                            pop_level_vec,
                            ncores=ncores)
    return(pvals)
end
################################################################################

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

function getReplicate(SLiM_file)
    fname = basename(SLiM_file)
    no_ext_fname = split(fname, ".")[1]
    replicate_field = split(no_ext_fname, "_")[5]
    replicate = split(replicate_field, "Rep")[2]
    return replicate
end


function GeneratePowerAnalysisDataFrame(ncores = 4)
    ## Add rows to this dataframe, one by one.
    power_analysis_df = DataFrame(MutatorStatus = String[],
                                  GeneModule = String[],
                                  Replicate = String[],
                                  Variable = String[],
                                  VariableValue = Float64[],
                                  Count = Float64[],
                                  Pval = Float64[])
    ## Variable parameters to test.
    timeseries_length_vec = [500, 1000, 2000, 5000]
    sampling_interval_vec = [25, 50, 100, 500]
    metagenomic_filtering_threshold_vec = [0, 0.001, 0.01, 0.05, 0.1]
    num_genes_from_module_vec = [1, 10, 25, 50, 100]
    
    outdir = "../results/SLiM-results"
    genome_metadata_csv_path = joinpath(outdir, "SLiM_geneIDs.csv")
    genome_metadata = CSV.read(genome_metadata_csv_path, DataFrame)
    
    gene_set_path_dict = Dict("positive" => joinpath(outdir, "SLiM_positive_module.csv"),
                              "neutral" => joinpath(outdir, "SLiM_neutral_module.csv"),
                              "purifying" => joinpath(outdir, "SLiM_purifying_module.csv"))
    gene_module_df_dict = Dict(k => CSV.read(v, DataFrame) for (k,v) in gene_set_path_dict)
    
    SLiM_dataset_dir = "../results/SLiM-results/SLiM-replicate-runs"
    ## skip over any irrelevant files (like .DS_Store on macs!)
    SLiM_dataset_vec = [x for x in readdir(SLiM_dataset_dir; join=true) if startswith(basename(x),"SLiM_")]
    
    for SLiM_file in SLiM_dataset_vec
        ## name as "Hypermutator" or "Nonmutator" based on the filename.
        pop_name = getMutationRate(SLiM_file)
        ## get the replicate from the filename.
        replicate = getReplicate(SLiM_file)
        
        for exp_length in timeseries_length_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; timeseries_length = exp_length)
            for (module_type, module_df) in gene_module_df_dict         
                pval_df = RunSTIMS_on_data(STIMS_input, genome_metadata, module_df, ncores)
                println(module_type, ":")
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, replicate,"timeseries_length", exp_length, pval_df.count[1], pval_df.pval[1]])
            end
        end
        
        for sampling_int in sampling_interval_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; sampling_interval = sampling_int)
            for (module_type, module_df) in gene_module_df_dict
                println(module_type, ":")
                pval_df = RunSTIMS_on_data(STIMS_input, genome_metadata, module_df, ncores)
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, replicate, "sampling_interval", sampling_int, pval_df.count[1], pval_df.pval[1]])
            end
        end
        
        for threshold in metagenomic_filtering_threshold_vec
            STIMS_input = SLiM_output_to_STIMS_input(
                SLiM_file, pop_name; freq_threshold = threshold)
            for (module_type, module_df) in gene_module_df_dict
                println(module_type, ":")
                pval_df = RunSTIMS_on_data(STIMS_input, genome_metadata, module_df, ncores)
                println(pval_df)
                push!(power_analysis_df, [pop_name,  module_type, replicate, "filtering_threshold", threshold, pval_df.count[1], pval_df.pval[1]])
            end
        end
        ## Note that this case is a bit different.
        for num_genes in num_genes_from_module_vec
            STIMS_input = SLiM_output_to_STIMS_input(SLiM_file, pop_name)
            for (module_type, module_df) in gene_module_df_dict
                ## Only run on the first num_genes in the module.
                filtered_gene_module_df = module_df[1:num_genes, :]
                println(module_type, ":")
                pval_df = RunSTIMS_on_data(STIMS_input, genome_metadata,
                                                 filtered_gene_module_df, ncores)
                println(pval_df)
                ## push a new row onto the power analysis DataFrame.
                push!(power_analysis_df, [pop_name,  module_type, replicate, "num_genes_from_module", num_genes, pval_df.count[1], pval_df.pval[1]])
            end
        end
    end
    return(power_analysis_df)
end


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "ncores"
        help = "number of cores for parallel computation."
        arg_type = Int64
        default = 4
    end

    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    ncores = parsed_args["ncores"]
    ## make the DataFrame for the power analysis
    power_analysis_df = GeneratePowerAnalysisDataFrame(ncores)
    ## and write to file. Make figures with analyze-STIMS-power.R.
    CSV.write("../results/SLiM-results/STIMS-power-analysis-data.csv", power_analysis_df)
end


## run the script.
main()
