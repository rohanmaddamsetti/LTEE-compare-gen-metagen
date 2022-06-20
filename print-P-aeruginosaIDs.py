#!/usr/bin/env python

'''
print-P-aeruginosaIDs.py by Rohan Maddamsetti

This particular reference genbank file is the one used by Mehta et al. (2018).

Usage: python print-P-aeruginosaIDs.py -i ../data/PAO11-AE004091.2.gbk > ../results/PAO11_IDs.csv

'''

import argparse
import re
from Bio import SeqIO
        
def main():
    parser = argparse.ArgumentParser(description='print csv file of CDS IDs in Genbank file.')

    parser.add_argument("-i", "--input", help='Input genbank file',required=True)
    
    args = parser.parse_args()

    ## now parse the genome for CDS.
    genome = next(SeqIO.parse(args.input, "genbank"))
    print(','.join(['Gene','PSEAE_locus_tag','gene_length','product', 'start', 'end', 'strand']))    
    for feat in genome.features:
        if feat.type != 'CDS': continue ## only consider protein-coding genes
        my_start = feat.location.start
        my_end = feat.location.end
        my_strand = feat.location.strand
        length = my_end - my_start
        locus_tag = feat.qualifiers['locus_tag'].pop()
        try:
            gene = feat.qualifiers['gene'].pop()
        except KeyError:
            gene = locus_tag
        try:
            product = feat.qualifiers['product'].pop()
        except:
            product = 'NA'
        ## strip all punctuation that could cause parsing problems.
        product = re.sub('[,;()]', '', product)
        ## start+1 (but end+0) to be consistent with genbank format.
        print(','.join([str(x) for x in (gene,locus_tag,length,product, my_start+1, my_end, my_strand)]))


main()
