import os
import re
import sys
import argparse
import time
from rnaseq_utils import RNAseqDataset
import bedgraph_utils as bgu

###################
# INPUT ARGUMENTS #
###################

parser = argparse.ArgumentParser()
parser.add_argument(
    "-N", "--name", dest='NAME',
    help="Prefix to add to the output filename(s)",
    default='bed12', type=str
)
parser.add_argument(
    "-D", "--directory", dest='DIRECTORY', type=str, default='.',
    help="Output directory for BEDGRAPH files."
)
parser.add_argument(
    "-O", "--output", dest='OUTPUT', type=str, default='stdout',
    help="Filepath to write multimapper-rescued BED12 file."
)
parser.add_argument(
    "-S","--stranded", dest='STRANDED',
    help="Filter the signal into multiple bedgraph files based on strand.",
    default=False, action='store_true'
)
parser.add_argument(
    "--sorted", dest='SORTED',
    help="If input BED file is sorted, output can be streamed.",
    default=False, action='store_true'
)
parser.add_argument(
    "--split_samples", dest='SPLIT_SAMPLES',
    help="Filter the signal into multiple columns in the bedgraph based on the sample name.",
    default=False, action='store_true'
)
parser.add_argument(
    "--multi_key", dest='MULTI_KEY',
    help="List of names that should be output as separate columns if --split_samples.",
    default=[], nargs='+'
)
parser.add_argument(
    "--label_includes", dest='LABEL_INCLUDES',
    help="Only output signal from lines that include this character string in the label.",
    default=None, type=str
)
parser.add_argument(
    "--type", dest='TYPE',
    help="Whether to output full read coverage, 5', or 3' end signal.",
    default='all', choices=['all','5','3']
)
parser.add_argument(
    "--bed", dest='BED',
    help="Input is in a 15-column BED format.",
    default=False, action="store_true"
)
parser.add_argument(
    "FILENAME", nargs='?'
)
args = parser.parse_args()

################################
# ENVIRONMENT SETUP: FUNCTIONS #
################################

def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]


def file_extension(filename):
    """Returns the file extension of a path as a string"""
    split_name = filename.split('.')
    if len(split_name) == 1:
        return None
    else:
        return split_name[-1].lower()

def input_is_valid(filename):
    """Boolean if the file's extension is valid (BED, ELR)"""
    if file_extension(filename) in ['bed','elr']:
        return True
    else:
        return False

def get_block_ranges(chromStart,blockStarts,blockSizes):
    """Reorganizes the BED12 block format into a list of doubles."""
    sizes = [int(i) for i in blockSizes.split(',')]
    starts = [int(i) for i in blockStarts.split(',')]
    s = int(chromStart)
    
    return [(s+start, s+start+size-1) for start,size in zip(starts,sizes)]

def output_lines(lines,output=args.OUTPUT):
    """Takes a list of bed lines and writes
    them to the output stream.
    """
    if output == 'stdout':
        for output_string in lines:
            print(output_string.rstrip())
    else:
        for output_string in lines:
            output.write('{}\n'.format(output_string.rstrip()))

def bed12_coverage_dict(bed_in,label_includes=None,stranded=False,coverage_type='all',split_samples=False):
    """Generates a bedgraph_dict from a BED12 file."""
    bedgraph_dict = {'+':{}, '-':{}, '.':{}}
    for line in bed_in:
        if line[0] == '#':
            continue
        
        chrom,start,end,readname,score,strand,mmnum,mmorder,rgb,blocknum,blocksizes,blockstarts,reads,samplename,label = line.rstrip().split('\t')
        # Check if the read should be filtered based on its label
        if label_includes is not None:
            if label_includes not in label:
                continue
        
        # Get the positions
        block_ranges = get_block_ranges(start,blockstarts,blocksizes)
        if coverage_type == 'all':
            positions_to_add = flatten([range(a,b+1) for a,b in block_ranges])
        elif coverage_type == '5':
            if strand == '+':
                positions_to_add = [block_ranges[0][0]]
            elif strand == '-':
                positions_to_add = [block_ranges[-1][-1]]
            else:
                positions_to_add = []
        elif coverage_type == '3':
            if strand == '+':
                positions_to_add = [block_ranges[-1][-1]]
            elif strand == '-':
                positions_to_add = [block_ranges[0][0]]
            else:
                positions_to_add = []
        
        if stranded:
            strand_to_add = strand
        else:
            strand_to_add = '.'
        
        if chrom not in bedgraph_dict[strand_to_add]:
            bedgraph_dict[strand_to_add][chrom] = {}
        
        for pos in positions_to_add:
            if split_samples:
                if samplename not in multi_key:
                    multi_key.append(samplename)
                if pos not in bedgraph_dict[strand_to_add][chrom]:
                    bedgraph_dict[strand_to_add][chrom][pos] = {}
                
                bedgraph_dict[strand_to_add][chrom][pos][samplename] = bedgraph_dict[strand_to_add][chrom][pos].get(samplename,float(0)) + float(reads)
            else:
                bedgraph_dict[strand_to_add][chrom][pos] = bedgraph_dict[strand_to_add][chrom].get(pos,float(0)) + float(reads)
    
    return bedgraph_dict

def vector_to_bedgraph_lines(value_vector, chrom, chromstart=0, digits=3, multi_key=[]):
    """Converts a vector of values to bedgraph format, yielding lines."""
    start = 0
    prevpos = 0
    prevcount = None
    position = None
    vpos = None
    # Iterate through every position with values in the dict
    
    for vpos in [k for k,v in enumerate(value_vector) if v != 0]:
        if multi_key:
            all_values = [str(round(value_vector[vpos].get(k,0),digits)) for k in multi_key]
            count = '\t'.join(all_values)
        else:
            count = round(value_vector[vpos],digits)
        
        if count != prevcount or int(vpos) > 1 + prevpos:
            # The newly encountered value is not a continuation
            # of the previous value. Write the old run and start another.
            if prevcount and prevcount != 0:
                line_to_write = '\t'.join(
                    [
                        str(i) for i in [
                            chrom,
                            start + chromstart,
                            prevpos + 1 + chromstart,
                            prevcount
                        ]
                    ]
                )
                
                yield line_to_write
            
            start = vpos
        
        prevcount = count
        prevpos = int(vpos)

    if vpos and prevcount and prevcount != 0:
        line_to_write = '\t'.join(
            [
                str(i) for i in [
                    chrom,
                    start + chromstart,
                    prevpos + 1 + chromstart,
                    prevcount
                ]
            ]
        )
        
        yield line_to_write

def dict_to_bedgraph_lines(values_dict, chrom, chromstart=0, digits=3, multi_key=[]):
    """Converts dict of values to bedgraph format"""
    start = 0
    prevpos = 0
    prevcount = None
    position = None
    # Iterate through every position with values in the dict
    for position in sorted(list(values_dict.keys())):
        if multi_key:
            all_values = [str(round(values_dict[position].get(k,0),digits)) for k in multi_key]
            count = '\t'.join(all_values)
        else:
            count = round(values_dict[position],digits)
        
        if count != prevcount or int(position) > 1 + prevpos:
            # The newly encountered value is not a continuation
            # of the previous value. Write the old run and start another.
            if prevcount and prevcount != 0:
                line_to_write = '\t'.join(
                    [
                        str(i) for i in [
                            chrom,
                            start + chromstart,
                            prevpos + chromstart + 1,
                            prevcount
                        ]
                    ]
                ) + '\n'
                
                yield line_to_write
            
            start = position
        prevcount = count
        prevpos = int(position)
    
    if position and prevcount and prevcount != 0:
        line_to_write = '\t'.join(
            [
                str(i) for i in [
                    chrom,
                    start + chromstart,
                    prevpos + chromstart + 1,
                    prevcount
                ]
            ]
        ) + '\n'
        
        yield line_to_write


####################
# PARSE INPUT FILE #
####################

s={'+':'plus','-':'minus','.':'ns'}

if args.FILENAME:
    if args.FILENAME.split('.')[-1].lower() not in ['bed','bed12','elr']:
        print("\nERROR: input file must be BED12 format.")
        parser.print_help()
        sys.exit(1)
    bed_in = open(args.FILENAME)
elif not sys.stdin.isatty():
    bed_in = sys.stdin
else:
    print("\nERROR: requires BED12 file as input.")
    parser.print_help()
    sys.exit(1)

output_files = {}
if args.OUTPUT != 'stdout':
    output_files['.'] = open(args.OUTPUT, 'w')
else:
    output_files['.'] = 'stdout'

if args.STRANDED:
    for strand in ['+','-']:
        filename = '{}/{}'.format(
            args.DIRECTORY,
            '.'.join([args.NAME,s[strand],'bedgraph'])
        )
        output_files[strand] = open(filename, 'w')
    
    if output_files['.'] == 'stdout':
        filename = '{}/{}'.format(
            args.DIRECTORY,
            '.'.join([args.NAME,s['.'],'bedgraph'])
        )
        output_files['.'] = open(filename, 'w')


# If --split_samples, multi_key will be used to store sample names
if args.MULTI_KEY:
    multi_key = args.MULTI_KEY
else:
    multi_key = []

if args.SORTED:
    bed12_coverage_stream(
        bed_in,
        label_includes=args.LABEL_INCLUDES,
        stranded=args.STRANDED,
        coverage_type=args.TYPE,
        split_samples=args.SPLIT_SAMPLES
    )
else:
    bedgraph_dict = bed12_coverage_dict(
        bed_in,
        label_includes=args.LABEL_INCLUDES,
        stranded=args.STRANDED,
        coverage_type=args.TYPE,
        split_samples=args.SPLIT_SAMPLES
    )
    for strand in ['+','-','.']:
        if bedgraph_dict[strand]:
            if args.STRANDED:
                output_filename = '{}/{}'.format(
                    args.DIRECTORY,
                    '.'.join([args.NAME,s[strand],'bedgraph'])
                )
            elif strand == '.':
                output_filename = '{}/{}'.format(
                    args.DIRECTORY,
                    '.'.join([args.NAME,'bedgraph'])
                )
            else:
                continue
            
            bgu.write_bedgraph_from_dict(bedgraph_dict[strand], output_filename, digits=3, multi_key=sorted(multi_key))




