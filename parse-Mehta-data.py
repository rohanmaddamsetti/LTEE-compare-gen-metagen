#!/usr/bin/env python

'''
parse-Mehta-data.py by Rohan Maddamsetti.

Prints the Mehta data in a nice format to analyze using STIMS.

Input is supplementary data from Mehta et al. (2018):
The Essential Role of Hypermutation in Rapid Adaptation to Antibiotic Stress

Usage: python parse-Mehta-data.py > ../results/Mehta2018-hypermutators.csv

'''

from os import listdir
from os.path import join, basename

output_header = "Population,MutationIndex,Position,Gene,t0"
print(output_header)

input_files = [x for x in listdir("../data") if "Mehta" in x and x.endswith('.csv')]
for f in input_files:
    sample = basename(f).split('.')[0].split('-')[1]
    with open(join("../data",f)) as fh:
        input_header_fields = {} ## dict of index to name of field.
        for i, line in enumerate(fh):
            line = line.strip()
            fields = line.split(',')
            if i == 0: ## populate input_header_fields and continue to the data
                input_header_fields = {k:x for k,x in enumerate(fields)}
                continue 
            mut_index = fields[0]
            position = fields[1]
            ## remove any square brackets in the gene field (these are deletions).
            gene = fields[-2].strip('[]')
            description = fields[-1]
            ## find the first day that the mutation is observed.
            percent_in_field = [('%' in x) for x in fields]
            try:
                idx_of_first_obs = percent_in_field.index(True)
            except:
                ## there is at least one row in the supplement
                ## that seems to be erroneously missing the first occurrence.
                ## ignore such cases.
                continue
            firstday = input_header_fields[idx_of_first_obs]
            assert firstday.startswith('day ')
            t0 = firstday.split()[-1]
            print(','.join([sample,mut_index,position,gene,t0]))
