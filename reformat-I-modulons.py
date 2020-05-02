#!/usr/bin/env python

'''
reformat-I-modulons.py by Rohan Maddamsetti

This script reads in the listing of genes in each I-modulon in
precise-db-repo/data/imodulon_gene_names.txt,
and reformats it into an R-friendly format for metagenomics-analysis.R.

Usage: python reformat-I-modulons.py > ../results/gene-modules/genes-to-I-modulons.csv
'''

print("Gene,I.modulon")
in_fh = open("../precise-db-repo/data/imodulon_gene_names.txt","r")
for i, line in enumerate(in_fh):
    if i == 0: continue
    line = line.strip()
    fields = line.split(',')
    Imodulon = fields[0]
    genes = fields[1:]
    for g in genes:
        print(','.join([g,Imodulon]))
