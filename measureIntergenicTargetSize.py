#!/usr/bin/env python

''' measureIntergenicTargetSize.py by Rohan Maddamsetti'''

intergenic_fasta_f = '../results/REL606_noncoding_seqs.fasta'

bp_count = 0
with open(intergenic_fasta_f) as fh:
    for l in fh:
        l = l.strip()
        if l.startswith('>'):
            continue
        else:
            bp_count = bp_count + len(l)

print('Length of intergenic regions:', bp_count)
