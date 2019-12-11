import numpy
import sys
from math import fabs
from random import sample

'''
measureTargetSize.py by Rohan Maddamsetti.

Substantial code taken from Benjamin Good's LTEE-metagenomic analysis code,
In particular from parse_file.py from that Github repo.

IMPORTANT TODO: gene numbers in this script and don't match exactly with
aerobic-anaerobic-metagenomics.R. 
I need to debug this for publication.

'''

base_table = {'A':'T','T':'A','G':'C','C':'G'}

codon_table = { 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R',
'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H', 'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!' }

# calculate number of synonymous opportunities for each codon
codon_synonymous_opportunity_table = {}
for codon in list(codon_table.keys()):
    codon_synonymous_opportunity_table[codon] = {}
    for i in range(0,3):
        codon_synonymous_opportunity_table[codon][i] = -1 # since G->G is by definition synonymous, but we don't want to count it
        codon_list = list(codon)
        for base in ['A','C','T','G']:
            codon_list[i]=base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_opportunity_table[codon][i]+=1

bases = set(['A','C','T','G'])
substitutions = []
for b1 in bases:
    for b2 in bases:
        if b2==b1:
            continue
        
        substitutions.append( '%s->%s' % (b1,b2) )
        
codon_synonymous_substitution_table = {}
codon_nonsynonymous_substitution_table = {}
for codon in list(codon_table.keys()):
    codon_synonymous_substitution_table[codon] = [[],[],[]]
    codon_nonsynonymous_substitution_table[codon] = [[],[],[]]
    
    for i in range(0,3):
        reference_base = codon[i]
        
        codon_list = list(codon)
        for derived_base in ['A','C','T','G']:
            if derived_base==reference_base:
                continue
            substitution = '%s->%s' % (reference_base, derived_base)
            codon_list[i]=derived_base
            new_codon = "".join(codon_list)
            if codon_table[codon]==codon_table[new_codon]:
                # synonymous!
                codon_synonymous_substitution_table[codon][i].append(substitution)
            else:
                codon_nonsynonymous_substitution_table[codon][i].append(substitution)


def parse_mask_list(filename="../data/REL606.L20.G15.P0.M35.mask.gd"):
    
    # Add masks calculated in Tenaillon et al (Nature, 2016)
    # Downloaded from barricklab/LTEE-Ecoli/reference/ github repository
    
    mask_start_positions = []
    mask_end_positions = []
    
    file = open(filename,"r")
    file.readline() # header
    for line in file:
        items = line.split()
        start = int(items[4])
        length = int(items[5])
        mask_start_positions.append(start)
        mask_end_positions.append(start+length-1)
    
    # Add masking of prophage elements (Methods section of Tenaillon et al (Nature, 2016))
    mask_start_positions.append(880528)
    mask_end_positions.append(904682)
    
    return numpy.array(mask_start_positions), numpy.array(mask_end_positions)

def is_repeat_masked(position, position_gene_map=None):

    if position_gene_map==None:
        position_gene_map,effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map()
        
    if position in position_gene_map:
        if position_gene_map[position]=='repeat':
            return True
            
    return False

def create_annotation_map(gene_data=None, repeat_data=None, mask_data=None):

    if gene_data==None:
        gene_data = parse_gene_list()
        repeat_data = parse_repeat_list()
        mask_data = parse_mask_list()
    

    locus_tags, gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    repeat_names, repeat_start_positions, repeat_end_positions, repeat_complements = repeat_data
    mask_start_positions, mask_end_positions = mask_data
    
    position_gene_map = {}
    gene_position_map = {}

    num_masked_sites = 0

    # first mark things that are repeats
    # this takes precedence over all other annotations
    for start,end in zip(repeat_start_positions,repeat_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position]='repeat'
                num_masked_sites+=1
    
    # then mark masked things
    for start,end in zip(mask_start_positions, mask_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position]='repeat'
                num_masked_sites+=1
    
    
    # then greedily annotate genes at remaining sites
    for gene_name,start,end in zip(gene_names,gene_start_positions,gene_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                position_gene_map[position] = gene_name
                if gene_name not in gene_position_map:
                    gene_position_map[gene_name]=[]
                gene_position_map[gene_name].append(position)
    
    # remove 'partial' genes that have < 10bp unmasked sites
    for gene_name in list(sorted(gene_position_map.keys())):
        if len(gene_position_map[gene_name]) < 10:
            for position in gene_position_map[gene_name]:
                position_gene_map[position] = 'repeat'
            del gene_position_map[gene_name]
    
    # count up number of synonymous opportunities
    effective_gene_synonymous_sites = {}
    effective_gene_nonsynonymous_sites = {}
    
    substitution_specific_synonymous_sites = {substitution: 0 for substitution in substitutions}
    substitution_specific_nonsynonymous_sites = {substitution: 0 for substitution in substitutions}

    for gene_name,start,end,gene_sequence,strand in zip(gene_names, gene_start_positions, gene_end_positions, gene_sequences, strands):
        
        if gene_name not in gene_position_map:
            continue
        
        if strand=='forward':
            oriented_gene_sequence = gene_sequence
        else:
            oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
        
        for position in gene_position_map[gene_name]:
        
            if gene_name not in effective_gene_synonymous_sites:
                effective_gene_synonymous_sites[gene_name]=0
                effective_gene_nonsynonymous_sites[gene_name]=0
        
            if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                pass
            
            else:
        
                # calculate position in gene
                if strand=='forward':
                    position_in_gene = position-start
                else:
                    position_in_gene = end-position
        
                # calculate codon start
                codon_start = int(position_in_gene/3)*3
                codon = gene_sequence[codon_start:codon_start+3] 
                position_in_codon = position_in_gene%3       
            
                #print gene_name, start, end, position, codon,position_in_codon
            
                effective_gene_synonymous_sites[gene_name] += codon_synonymous_opportunity_table[codon][position_in_codon]/3.0
                effective_gene_nonsynonymous_sites[gene_name] += 1-codon_synonymous_opportunity_table[codon][position_in_codon]/3.0 
                
                for substitution in codon_synonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_synonymous_sites[substitution] += 1
                    
                for substitution in codon_nonsynonymous_substitution_table[codon][position_in_codon]:
                    substitution_specific_nonsynonymous_sites[substitution] += 1

    substitution_specific_synonymous_fraction = {substitution: substitution_specific_synonymous_sites[substitution]*1.0/(substitution_specific_synonymous_sites[substitution]+substitution_specific_nonsynonymous_sites[substitution]) for substitution in list(substitution_specific_synonymous_sites.keys())}
                    
    # then annotate promoter regions at remaining sites
    for gene_name,start,end in zip(gene_names,promoter_start_positions,promoter_end_positions):
        for position in range(start,end+1):
            if position not in position_gene_map:
                # position hasn't been annotated yet
                
                if gene_name not in gene_position_map:
                    # the gene itself has not been annotated
                    # so don't annotate the promoter
                    pass
                else:
                    position_gene_map[position] = gene_name
                    gene_position_map[gene_name].append(position)
    
    # calculate effective gene lengths
    effective_gene_lengths = {gene_name: len(gene_position_map[gene_name])-effective_gene_synonymous_sites[gene_name] for gene_name in list(gene_position_map.keys())}
    effective_gene_lengths['synonymous'] = sum([effective_gene_synonymous_sites[gene_name] for gene_name in list(gene_position_map.keys())])
    effective_gene_lengths['nonsynonymous'] = sum([effective_gene_nonsynonymous_sites[gene_name] for gene_name in list(gene_position_map.keys())])
    effective_gene_lengths['masked'] = num_masked_sites
    effective_gene_lengths['noncoding'] = calculate_genome_length()-effective_gene_lengths['synonymous']-effective_gene_lengths['nonsynonymous']-effective_gene_lengths['masked']
    
    return position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction      

def calculate_synonymous_nonsynonymous_target_sizes():
    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction  = create_annotation_map()
    return effective_gene_lengths['synonymous'], effective_gene_lengths['nonsynonymous'], substitution_specific_synonymous_fraction   

def calculate_reverse_complement_sequence(dna_sequence):
    return "".join(base_table[base] for base in dna_sequence[::-1])
 
def calculate_codon_sequence(dna_sequence):
    return "".join(codon_table[dna_sequence[3*i:3*i+3]] for i in range(0,len(dna_sequence)/3))

def get_closest_gene_name(location, gene_data):
    gene_names, start_positions, end_positions, gene_sequences, strands = gene_data

    closest_start_idx = numpy.fabs(start_positions-location).argmin()
    closest_end_idx = numpy.fabs(end_positions-location).argmin()
    if fabs(start_positions[closest_start_idx]-location) < fabs(end_positions[closest_end_idx]-location):
        return gene_names[closest_start_idx]
    else:
        return gene_names[closest_end_idx]

var_types = ['synonymous','missense','nonsense','noncoding','indel','sv']

def annotate_gene(position, position_gene_map):
    
    if position in position_gene_map:
        gene_name = position_gene_map[position]
    else:
        gene_name = 'intergenic'
        
    return gene_name

def annotate_variant(position, allele, gene_data, position_gene_map):
    
    gene_names, gene_start_positions, gene_end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
    # get gene
    gene_name = annotate_gene(position, position_gene_map)
    
    if allele.startswith('Depth'):
        var_type = 'unknown'
    elif allele.startswith('MOB') or allele.startswith('junction'): 
        var_type = 'sv'
    elif allele.startswith('indel'):
            var_type = 'indel'
    elif allele[1:3]=='->':
        # a SNP, so annotate it
        if gene_name=='intergenic':
            var_type = 'noncoding'
        elif gene_name=='repeat':
            var_type = 'repeat'
        else:
            # must be in a real gene
            # so get it
            i = gene_names.index(gene_name)
            
            gene_start_position = gene_start_positions[i]
            gene_end_position = gene_end_positions[i]
            promoter_start_position = promoter_start_positions[i]
            promoter_end_position = promoter_end_positions[i]
            gene_sequence = gene_sequences[i]
            strand = strands[i]
            
            if position<gene_start_position or position>gene_end_position:
                #var_type='promoter'
                var_type='noncoding' # (promoter)
            else:
            
                if gene_name.startswith('tRNA') or gene_name.startswith('rRNA'):
                    var_type='noncoding'
                else:
                
                    # calculate position in gene
                    if strand=='forward':
                        position_in_gene = position-gene_start_position
                        oriented_gene_sequence = gene_sequence
                        new_base = allele[3]
                    else:
                        position_in_gene = gene_end_position-position
                        oriented_gene_sequence = calculate_reverse_complement_sequence(gene_sequence)
                        new_base = base_table[allele[3]]
            
                
                    # calculate codon start
                    codon_start = int(position_in_gene/3)*3
                    codon = oriented_gene_sequence[codon_start:codon_start+3]
                    codon_list = list(codon) 
                    position_in_codon = position_in_gene%3 
                    codon_list[position_in_codon]=new_base
                    new_codon="".join(codon_list)
                    if codon_table[codon]==codon_table[new_codon]:
                        var_type='synonymous'
                    else:
            
                        if codon_table[new_codon]=='!':
                            var_type='nonsense'
                        else:
                            var_type='missense'
    else:
        sys.stderr.write("Unknown: %s\n" % allele)
        var_type='unknown'
    
    return gene_name, var_type
        

def get_closest_repeat_idx(location, repeat_data):
    repeat_names, start_positions, end_positions, complements = repeat_data

    closest_start_idx = numpy.fabs(start_positions-location).argmin()
    closest_end_idx = numpy.fabs(end_positions-location).argmin()
    if fabs(start_positions[closest_start_idx]-location) < fabs(end_positions[closest_end_idx]-location):
        return closest_start_idx
    else:
        return closest_end_idx

def get_repeat_idx(location, repeat_data):
    repeat_names, start_positions, end_positions, complements = repeat_data
    
    repeat_idxs = numpy.nonzero((location >= start_positions-100) * (location <= end_positions+100))[0]
    if len(repeat_idxs) > 0:
        return repeat_idxs[0]
    else:
        return -1
        

def in_repeat_region(location, repeat_data):
    
    repeat_names, start_positions, end_positions, complements = repeat_data
    
    repeat_idxs = numpy.nonzero((location >= start_positions) * (location <= end_positions))[0]
    if len(repeat_idxs) > 0:
        return True
    else:
        return False
        

def parse_reference_genome(filename="../data/REL606.7.gbk"):
    reference_sequences = []
    
    # GBK file
    if filename[-3:] == 'gbk':
        file = open(filename,"r")
        origin_reached = False
        for line in file:
            if line.startswith("ORIGIN"):
                origin_reached=True
            if origin_reached:
                items = line.split()
                if items[0].isdigit():
                    reference_sequences.extend(items[1:])    
        file.close()
    
    # FASTA file
    else:
        file = open(filename,"r")
        file.readline() # header
        for line in file:
            reference_sequences.append(line.strip())
        file.close()
    
    reference_sequence = "".join(reference_sequences).upper()
    return reference_sequence

def calculate_genome_length(reference_sequence=None):
    if reference_sequence==None:
        reference_sequence=parse_reference_genome()
    return len(reference_sequence)
    
def print_reference_fasta(reference_sequence):
    print(">chr1")
    for i in range(0,len(reference_sequence),70):
        print(reference_sequence[i:min([len(reference_sequence),i+70])])

def create_gene_size_map(effective_gene_lengths=None):

    if effective_gene_lengths==None:
        reference_sequence = parse_reference_genome()
        gene_data = parse_gene_list()
        repeat_data = parse_repeat_list()
        mask_data = parse_mask_list()
        gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    
        position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(gene_data, repeat_data, mask_data)

    
    excluded_genes=set(['synonymous','nonsynonymous','noncoding','masked'])
    
    gene_size_map = {}
    for gene_name in list(effective_gene_lengths.keys()):
        
        #if gene_name.startswith('tRNA'):
        #    print gene_name 
            
        if gene_name in excluded_genes:
            continue
            
        gene_size_map[gene_name] = effective_gene_lengths[gene_name]
        
    return gene_size_map    

#####################################################################
#
# Reads through the Genbank file for the reference and 
# compiles a list of genes, tRNAs, etc.
#
#####################################################################        
def parse_gene_list(reference_sequence=None, filename="../data/REL606.7.gbk"):

    features = set(['CDS','gene','tRNA','rRNA','repeat_region'])

    if reference_sequence==None:
        reference_sequence=parse_reference_genome()

    observed_gene_names = set()

    locus_tags = []
    gene_names = []
    start_positions = []
    end_positions = []
    promoter_start_positions = []
    promoter_end_positions = []
    gene_sequences = []
    strands = []

    file = open(filename,"r")
    line = file.readline()
    while line!="":
        
        items = line.split()
        feature = items[0]

        locus_name = ""
        gene_name = ""
        feature_location = ""
        
        if feature=='CDS':
            feature_location = items[1]
        
            line = file.readline().strip()
            
            gene_name=""
            locus_name=""
            is_pseudo=False
            
            while line.split()[0] not in features:
                
                if line.startswith('/gene'):
                    gene_name = line.split('=')[1].strip()[1:-1]
                if line.startswith('/locus_tag'):
                    locus_name = line.split('=')[1].strip()[1:-1] 
                if line.startswith('/pseudo'):
                    is_pseudo=True   
        
                line = file.readline().strip()
                
            if gene_name=="":
                gene_name = locus_name
                
            if is_pseudo:
                gene_name = ""
                     
            # done here
          
        elif feature=='tRNA' or feature=='rRNA':
        
            feature_location = items[1]
            #print feature_location
        
            while not line.strip().startswith('/gene'):
                line = file.readline().strip()
            gene_name = line.split('=')[1].strip()[1:-1]
            gene_name = '%s:%s' % (feature, gene_name)
            
            
        else:
            # nothing to see here
            line = file.readline().strip()
        
        # If the element has a feature location string and a name
        # it should either be a gene, tRNA, or rRNA, so let's get details
        if feature_location!="" and gene_name!="":
        
            location_str = feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
            location_strs = location_str.split(",")
            
            for location_str in location_strs:
            
                locations = [int(subitem) for subitem in location_str.split("..")]
            
                gene_start = locations[0]
                gene_end = locations[1]
                
                if feature=="CDS":
                    gene_sequence = reference_sequence[gene_start-1:gene_end]
                else:
                    gene_sequence = ""
                
                strand = 'forward'
                promoter_start = gene_start - 100 # by arbitrary definition, we treat the 100bp upstream as promoters
                promoter_end = gene_start - 1
                        
                
                if gene_sequence!="" and (not len(gene_sequence)%3==0):
                    pass
                    print(gene_name, gene_start, "Not a multiple of 3")
                    
                if feature_location.startswith('complement'):
                    strand='reverse'
                    promoter_start = gene_end+1
                    promoter_end = gene_end+100
                    
                    if gene_sequence=="":
                        promoter_end = promoter_start
                
                
                # record information
                
                # first make sure gene name is unique
                i = 1
                old_gene_name = gene_name
                while gene_name in observed_gene_names:
                    i+=1
                    gene_name = "%s_%d" % (old_gene_name,i)

                ## make sure that the locus_tag is defined,
                ## and that it is not already in the list.
                if not len(locus_name):
                    continue
                if locus_name in locus_tags:
                    continue
                    
                locus_tags.append(locus_name)
                start_positions.append(gene_start)
                end_positions.append(gene_end)
                promoter_start_positions.append(promoter_start)
                promoter_end_positions.append(promoter_end)
                gene_names.append(gene_name)
                gene_sequences.append(gene_sequence)
                strands.append(strand)
                observed_gene_names.add(gene_name)
        
    file.close()
    
    # sort genes based on start position
    
    locus_tags, gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = (list(x) for x in zip(*sorted(zip(locus_tags, gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands), key=lambda pair: pair[2])))
    
    return locus_tags, gene_names, numpy.array(start_positions), numpy.array(end_positions), numpy.array(promoter_start_positions), numpy.array(promoter_end_positions), gene_sequences, strands
    
def parse_repeat_list(filename="../data/REL606.7.gbk"):

    repeat_names = []
    start_positions = []
    end_positions = []
    complements = []
 
    file = open(filename,"r")
    line = file.readline()
    while line!="":
        items = line.split()
        feature = items[0]
        
        if feature == 'repeat_region':
            feature_location = items[1]
        
            # Get name of mobile element
            repeat_name = 'unknown'
            
            line = file.readline()
            while line.strip()[0]=='/':
                if line.strip().startswith('/mobile_element'):
                    repeat_name = line.split('=')[1].strip()[1:-1]
                line = file.readline()
            
            # Finished at next non '/' entry, presumably next feature

            if feature_location.startswith('complement'):
                complement = True
            else:
                complement = False            
            
            location_str = feature_location.lstrip("complement(").lstrip("join(").rstrip(")")
            location_strs = location_str.split(",")
            for location_str in location_strs:
            
                locations = [int(subitem) for subitem in location_str.split("..")]
                start_positions.append(locations[0])
                end_positions.append(locations[1])
                repeat_names.append(repeat_name)
                complements.append(complement)
        
        else:
            line = file.readline()
    file.close()
    
    return repeat_names, numpy.array(start_positions), numpy.array(end_positions), complements
        
def count_differences(mutation_list_1, mutation_list_2):
    unique_mutations = set()
    unique_mutations.update(mutation_list_1)
    unique_mutations.update(mutation_list_2)
    return 2*len(unique_mutations)-len(mutation_list_1)-len(mutation_list_2)

def parse_specific_tags(gene_f):
    tags = []
    tags_fh = open(gene_f,'r')
    for i,line in enumerate(tags_fh):
        line = line.strip()
        if i == 0: ## skip header.
            continue
        fields = line.split(',')
        tag = fields[3].strip("\"")
        ## don't append duplicates or empty entries.
        if len(tag) and tag not in tags:
            tags.append(tag)
    return tags

def parse_anaerobic_tags():
    return parse_specific_tags("../results/anaerobic-specific-genes.csv")

def parse_aerobic_tags():
    return parse_specific_tags("../results/aerobic-specific-genes.csv")

def filter_gene_data(gene_data, tags):
    locus_tags, gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data

    indices = [i for i in range(len(locus_tags)) if locus_tags[i] in tags]
    filtered_locus_tags = [locus_tags[i] for i in indices]
    filtered_gene_names = [gene_names[i] for i in indices]
    filtered_start_positions = [start_positions[i] for i in indices]
    filtered_end_positions = [end_positions[i] for i in indices]
    filtered_promoter_start_positions = [promoter_start_positions[i] for i in indices]
    filtered_promoter_end_positions = [promoter_end_positions[i] for i in indices]
    filtered_gene_sequences = [gene_sequences[i] for i in indices]
    filtered_strands = [strands[i] for i in indices]

    filtered_gene_data = (filtered_locus_tags, filtered_gene_names, filtered_start_positions, filtered_end_positions, filtered_promoter_start_positions, filtered_promoter_end_positions, filtered_gene_sequences, filtered_strands)
    
    return(filtered_gene_data)

def filter_anaerobic_gene_data(gene_data):
    return filter_gene_data(gene_data, parse_anaerobic_tags())

def filter_aerobic_gene_data(gene_data):
    return filter_gene_data(gene_data, parse_aerobic_tags())

def choose_random_tags(gene_data,N):
    ''' choose N genes out of gene_data.
    return: a list of randomly selected locus_tags.'''
    locus_tags = gene_data[0]
    random_tags = sample(locus_tags, N)
    return(random_tags)

def tag_to_gene(tags,gene_data):
    locus_tags, gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = gene_data
    tag_dict = {}
    for i,l in enumerate(locus_tags):
        if l in tags:
            tag_dict[l] = gene_names[i]
    mapped_genes = [tag_dict[x] for x in tags]
    return mapped_genes

def print_random_gene_csv(random_tags,gene_data, outf):
    random_genes = tag_to_gene(random_tags, gene_data)
    outfh = open(outf,'w')
    outfh.write('gene,REL606_locus_tag\n')
    for g,t in zip(random_genes,random_tags):
        outfh.write(g+','+t+'\n')
    outfh.close()
    
def print_target_size_statistics(setname, data):
    ''' 
    to write out a table of the results I need in
    aerobic-anaerobic-metagenomics.R, represent rows
    in the csv file as a list of strings.

    Header: 
    set, size, total_gene_length, synon_sites, non_synon_sites
    Returns a list of strings.
    '''

    position_gene_map, effective_gene_lengths, substitution_specific_synonymous_fraction = create_annotation_map(data, repeat_data, mask_data)

    locus_tags, gene_names, start_positions, end_positions, promoter_start_positions, promoter_end_positions, gene_sequences, strands = data
    total_gene_set_length = sum([fabs(end_pos-start_pos) for start_pos, end_pos in zip(start_positions,end_positions)])
    
    print("Masked:", effective_gene_lengths['masked'])
    print("Synonymous sites:", effective_gene_lengths['synonymous'])
    print("Nonsynonymous sites:", effective_gene_lengths['nonsynonymous'])
    print("Noncoding sites:", effective_gene_lengths['noncoding'])
    
    print("Nonsynonymous:synonymous ratio:", effective_gene_lengths['nonsynonymous']/effective_gene_lengths['synonymous'])
    print("Noncoding:synonymous ratio:", effective_gene_lengths['noncoding']/effective_gene_lengths['synonymous'])

    print('Total number of sites:', total_gene_set_length)
    print(len(data[0]), "genes")

    set_size = len(data[0])
    synon_sites = effective_gene_lengths['synonymous']
    nonsynon_sites = effective_gene_lengths['nonsynonymous']
    
    fields = [str(x) for x in (setname,set_size,total_gene_set_length,synon_sites,nonsynon_sites)]
    return ','.join(fields)

def write_stats_to_file(f,lines):
    header = 'set, size, total_gene_length, synon_sites, non_synon_sites'
    with open(f,'w') as fh:
        fh.write(header+'\n')
        for l in lines:
            fh.write(l+'\n')
    fh.close()
    
if __name__=='__main__':
    reference_sequence = parse_reference_genome()
    gene_data = parse_gene_list()
    repeat_data = parse_repeat_list()
    mask_data = parse_mask_list()

    ## filter gene_data based on whether in anaerobic-specific genes.
    anaerobic_gene_data = filter_anaerobic_gene_data(gene_data)
    print('\nANAEROBIC-SPECIFIC GENE STATISTICS:\n')
    anaerobic_string = print_target_size_statistics('anaerobic', anaerobic_gene_data)
    
    ## filter gene_data based on whether in aerobic-specific genes.
    aerobic_gene_data = filter_aerobic_gene_data(gene_data)
    print('\nAEROBIC-SPECIFIC GENE STATISTICS:\n')
    aerobic_string = print_target_size_statistics('aerobic', aerobic_gene_data)

    ## now print statistics for the whole genome.
    print('\nWHOLE GENOME STATISTICS:\n')
    print("Total:", len(reference_sequence))
    total_genome_string = print_target_size_statistics('genome', gene_data)

    print('\nSTATISTICS FOR A SET OF RANDOM GENES, THE SIZE OF ANAEROBIC GENE SET: ')
    ## now print statistics for a random set of genes with the same length as the anaerobic genes.

    anaerobic_rando_set = choose_random_tags(gene_data,len(anaerobic_gene_data[0]))
    anaerobic_random_csvfile = '../results/random-anaerobic-set.csv'
    print('printing random genes to:',anaerobic_random_csvfile)
    print_random_gene_csv(anaerobic_rando_set,gene_data, anaerobic_random_csvfile)
    anaerobic_rando_gene_data = filter_gene_data(gene_data, anaerobic_rando_set)
    random_anaerobic_string = print_target_size_statistics('random_anaerobic', anaerobic_rando_gene_data)

    print('\nSTATISTICS FOR A SET OF RANDOM GENES, THE SIZE OF AEROBIC GENE SET: ')
    ## now print statistics for a random set of genes with the same length as the anaerobic genes.
    aerobic_rando_set = choose_random_tags(gene_data,len(aerobic_gene_data[0]))
    aerobic_random_csvfile = '../results/random-aerobic-set.csv'
    print('printing random genes to:',aerobic_random_csvfile)
    print_random_gene_csv(aerobic_rando_set,gene_data, aerobic_random_csvfile)
    aerobic_rando_gene_data = filter_gene_data(gene_data,aerobic_rando_set)
    random_aerobic_string = print_target_size_statistics('random_aerobic', aerobic_rando_gene_data)    

    strings_to_write = [anaerobic_string, aerobic_string, random_anaerobic_string, random_aerobic_string, total_genome_string]

    statoutf = '../results/target_size.csv'
    print('writing statistics for aerobic-anaerobic-metagenomics.R to:',statoutf)
    write_stats_to_file(statoutf, strings_to_write)
    
