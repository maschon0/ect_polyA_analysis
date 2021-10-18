import sys
import os
import re
import math
import copy
import argparse
import fasta_utils as fu
import rnaseq_utils as ru
import numpy as np
from collections import Counter

desc = (
    "Calculates the BEDGRAPH signal overlapping a BED/BED12 feature."
)
parser = argparse.ArgumentParser(description=desc)

# add arguments to the ArgumentParser
parser.add_argument(
    '-G', '--genome', dest='genome', type=str, 
    help='input genome fasta file',
    required=True
)
parser.add_argument(
    '--features', dest='features', type=str,
    help='3-prime feature annotation file',
    required=True
)
parser.add_argument(
    '--extend', dest='extend', type=int,
    help='Number of nucletides to extend 3-prime end of annotation',
    default=0
)
parser.add_argument(
    '-I', '--input', dest='input', type=str,
    help='input bedgraph coverage files (+,-)',
    required=True, nargs=2
)
parser.add_argument(
    '--type', dest='type', type=str,
    help='', default='consensus',
    choices=['consensus','median']
)


args = parser.parse_args()

def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def get_consensus(numpy_array, ranges, strand):
    """Returns the (absolute) position in a numpy array
    that satisfies two rules within the set of (start, end) ranges:
      1) The position with the most reads
      2) The more downstream (strand-dependent) position in ties
    """
    sub_array = np.concatenate([numpy_array[start:end] for start,end in ranges])
    traceback = np.concatenate([range(start,end) for start,end in ranges])  
    total = sum([int(i) for i in sub_array])
    if total == 0:
        return None, 0
    
    hits = np.where(sub_array == np.max(sub_array))[0]
    if strand == '+':
        return traceback[max(hits)], total
    else:
        return traceback[min(hits)], total

def get_median(numpy_array, ranges, strand):
    """Returns the median (absolute) position in a set of ranges within an array.
    Odd-numbered positions are rounded downstream based on strand.
    """
    sub_array = np.concatenate([numpy_array[start:end] for start,end in ranges])
    traceback = np.concatenate([range(start,end) for start,end in ranges])
    total = sum([int(i) for i in sub_array])
    if total == 0:
        return None, 0
    
    positional_hits = flatten([[pos]*int(i) for pos,i in enumerate(sub_array)])
    if len(positional_hits) == 1:
        middle = 0
    else:
        if strand == '+':
            middle = math.floor(len(positional_hits)/2)
        else:
            middle = math.floor((len(positional_hits)-1)/2)
    
    return traceback[positional_hits[middle]], total

if __name__=='__main__':
    genome = fu.import_genome(args.genome)
    chromosomes = dict([(k,len(v)) for k,v in genome.items()])
    
    coverage = {}
    ingraphs = args.input
    assert len(ingraphs) == 2
    for graph in ingraphs:
        if graph == ingraphs[0]:
            strand = '+'
        else:
            strand = '-'
        
        coverage[strand] = dict([(k,np.zeros(v,dtype=np.float32)) for k,v in chromosomes.items()])
        coverage_file = open(graph)
        for line in coverage_file:
            chrom,start,end,count = line.rstrip().split()
            if chrom in coverage[strand]:
                coverage[strand][chrom][int(start):int(end)] += float(count)
    
    bedfile=open(args.features)
    BODY = ru.RNAseqDataset(chrom_array=sorted(list(chromosomes.keys())),record_names=True)
    for line in bedfile:
        if line[0]=='#':
            continue
        
        # Add this exon to the dict
        BODY.add_read_from_BED(line,gaps_are_junctions=False)
    
    bedfile.close()
    position = None
    for name,read in BODY.read_names.items():
        # Get the dominant feature
        chrom = BODY.chrom_array[read.chrom]
        strand = read.strand
        if strand == '+':
            ranges = [(a,b+args.extend) for a,b in read.ranges]
        elif strand == '-':
            ranges = [(a-args.extend,b) for a,b in read.ranges]
        
        if args.type == 'consensus':
            position, total = get_consensus(coverage[strand][chrom], ranges, strand)
        elif args.type == 'median':
            position, total = get_median(coverage[strand][chrom], ranges, strand)
        
        
        if position is not None:
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(
                chrom, position, position+1, name, round(total, 1), strand
            ))


