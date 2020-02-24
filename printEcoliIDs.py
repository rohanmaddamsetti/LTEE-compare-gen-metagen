#!/usr/bin/env python

'''
Usage: python printEcoliIDs.py -i ../data/REL606.7.gbk --mask ../data/REL606.L20.G15.P0.M35.RM-edited.mask.gd > ../results/REL606_IDs.csv

This script skips over repetitive regions, per Section 4.3.1
"Removing mutations in repetitive regions of the genome"
in the Supplemental Information of Ben Good's LTEE metagenomics paper.

We excluded mutations that arose in repetitive regions of the genome, as these can be difficult to resolve using short-read data. A site was marked as repetitive if (1) it was annotated as a repeat region in the REL606 reference, (2) it was present in the set of masked regions compiled by Tenaillon et al. [3], or (3) it fell within the putative prophage element identified by Tenaillon et al. [3] (REL606 genome coordinates 880528–904682).

I added this last region to Tenaillon's mask file, and saved it as:
REL606.L20.G15.P0.M35.RM-edited.mask.gd.

THE USER MUST MAKE SURE THAT THE MASK FILE CORRESPONDS TO THE GENBANK FILE!!!
'''

import argparse
import re
from Bio import SeqIO

def overlaps_any_region(cds_start, cds_end, regions):
    for r in regions:
        if overlaps_this_region(cds_start, cds_end, r):
            return True
    return False
    
def overlaps_this_region(cds_start, cds_end, region):
    ''' returns True if the intervals spanned by the region and the cds 
        overlap at all. returns False if there's no overlap. '''
    region_start, region_end = region
    assert region_start < region_end
    assert cds_start < cds_end
    if (cds_end < region_start) or (cds_start > region_end):
        return False
    return True

def is_in_any_masked_region(cds_start, cds_end, masked_regions):
    for mask_start, mask_len in masked_regions:
        mask_end = mask_start + mask_len
        if overlaps_this_region(cds_start, cds_end, (mask_start, mask_end)):
            return True
    return False

def get_masked_regions(maskfile):
    ''' 
    input: a mask file.
    output: a list of tuples (region_start, region_end).
    '''
    mask_f = open(maskfile, "r")
    masked_regions = []
    for l in mask_f:
        if l.startswith('#'): continue
        assert l.startswith('MASK')
        l = l.strip()
        ldata = l.split()
        ## Genbank starts from 1, Biopython starts from 0.
        region_start = int(ldata[4]) - 1 ## so subtract 1 from the start.
        region_len = int(ldata[5])
        region_tuple = (region_start, region_len)
        masked_regions.append(region_tuple)
    return masked_regions

def get_repeat_regions(genbankf):
    genome = next(SeqIO.parse(genbankf, "genbank"))
    repeats = []
    for feat in genome.features:
        ## misc_features are prophage that Jeff Barrick annotated.
        ## let's keep misc_features in, for now. But consider skipping.
        ##if (feat.type == "repeat_region") or (feat.type == "misc_feature"):
        if (feat.type == "repeat_region"):
            repeat_start = feat.location.start
            repeat_end = feat.location.end
            repeat_tuple = (repeat_start, repeat_end)
            repeats.append(repeat_tuple)
    return repeats
        

def main():
    parser = argparse.ArgumentParser(description='print csv file of CDS IDs in Genbank file, as long as they do not overlap a repetitive region or prophage.')

    parser.add_argument("-i", "--input", help='Input genbank file',required=True)
    parser.add_argument("-m", "--mask", help="genomediff mask file", required=False)
    
    args = parser.parse_args()

    masked_regions = () ## initialize to an empty tuple.
    if args.mask: ## then parse the mask file.
        masked_regions = get_masked_regions(args.mask)
    ## parse the genome for repetitive regions.
    annotated_repeats = get_repeat_regions(args.input)

    ## now parse the genome for CDS.
    genome = next(SeqIO.parse(args.input, "genbank"))
    print(','.join(['Gene','locus_tag','blattner','gene_length','product', 'start', 'end', 'strand']))    
    for feat in genome.features:
        if feat.type != 'CDS': continue ## only consider protein-coding genes
        my_start = feat.location.start
        my_end = feat.location.end
        my_strand = feat.location.strand
        if overlaps_any_region(my_start, my_end, annotated_repeats): continue
        if is_in_any_masked_region(my_start, my_end, masked_regions): continue
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
        try:
            note = feat.qualifiers['note'].pop()
            blattner = note.split('MG1655 equivalent: ')[-1].split()[-1]
            blattner = re.sub('[,;()]', '', blattner)
        except:
            blattner = 'NA'

        ## strip all punctuation that could cause parsing problems.
        product = re.sub('[,;()]', '', product)
        ## I should not have to do this next-- fix these corner cases later.
        blattner = re.sub('[,;()]', '', blattner)
        ## start+1 (but end+0) to be consistent with genbank format.
        print(','.join([str(x) for x in (gene,locus_tag,blattner,length,product, my_start+1, my_end, my_strand)]))


main()
