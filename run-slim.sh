#!/bin/bash

## run-slim.sh by Rohan Maddamsetti
## run SLiM at Ne = 1e7 on DCC.

sbatch -t 96:00:00 --mem=12G --wrap="slim WF-evolution-model-v12.slim"
