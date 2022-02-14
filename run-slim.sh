#!/bin/bash

## run-slim.sh by Rohan Maddamsetti
## run replicate non-mutator and hypermutator SLiM simulations on the HPC.

## non-mutator runs, for 5000 generations.
sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep1.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep2.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep3.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep4.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep5.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep6.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep7.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep8.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep9.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000_Rep10.txt\' bacterial-WF-model.slim"

## hypermutator runs, for 5000 generations.
sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep1.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep2.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep3.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep4.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep5.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep6.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep7.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep8.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep9.txt\' bacterial-WF-model.slim"

sbatch -t 8:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000_Rep10.txt\' bacterial-WF-model.slim"
