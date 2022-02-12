#!/bin/bash

## run-slim.sh by Rohan Maddamsetti
## run SLiM with varying parameters on HPC.

## non-mutator run, for 5000 generations.
sbatch -t 24:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-10 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-10_numgens5000.txt\' bacterial-WF-model.slim"

## hypermutator run, for 5000 generations.
sbatch -t 24:00:00 --mem=12G --wrap="slim -l -d Ne=1e6 -d mu=1e-8 -d numgens=5000 -d outfile=\'SLiM_Ne1000000_mu10-8_numgens5000.txt\' bacterial-WF-model.slim"
