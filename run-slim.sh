#!/bin/bash

## run-slim.sh by Rohan Maddamsetti
## run SLiM with varying parameters on HPC.

sbatch -t 96:00:00 --mem=12G --wrap="slim -l -d Ne=1e5 -d mu=1e-8 -d numgens=1000 -d outfile=\'SLiM_Ne10000_mu10-8_numgens1000.txt\' bacterial-WF-model.slim"
