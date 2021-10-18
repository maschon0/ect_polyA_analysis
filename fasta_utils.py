import re

IUPAChash = {}
IUPACnum = {}
IUPACbin = bytearray(b'\0')*256
IUPACcomp = {}
IUPACdeg = {}
# IUPAC code for nucleotides and ambiguities
keys = [
    'A','C','G','T','U',
    'M','R','W','S','Y','K',
    'V','H','D','B','N','X'
]
# REGEX representation of IUPAC nucleotides
values = [
    'A','C','G','T','T',
    '[AC]','[AG]','[AT]','[CG]','[CT]','[GT]',
    '[ACG]','[ACT]','[AGT]','[CGT]','.','-'
]

# integer values for a base-2 translation of A(0001), C(0010), G(0100), T(1000)
A = 1
C = 2
G = 4
T = 8
binaries = [
    A,C,G,T,T,
    A|C,A|G,A|T,C|G,C|T,G|T,
    A|C|G,A|C|T,A|G|T,C|G|T,A|C|G|T,0
]

# An ordered index of IUPAC characters by their binary value
bin_to_IUPAC = [
    'X','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'
]

# IUPAC representation of the complement of 'keys'
complements = [
    'T','G','C','A','A',
    'K','R','W','S','Y','M',
    'B','D','H','V','N','X'
]
# Number of possible nucleotides represented by letter
degeneracy = [
    1,1,1,1,1,
    2,2,2,2,2,2,
    3,3,3,3,4,0
]

for k,v,b,c,d in zip(keys, values, binaries, complements, degeneracy):
     IUPAChash[k] = v
     IUPACnum[k] = b
     IUPACbin[ord(k)] = b
     IUPACbin[ord(k.lower())] = b
     IUPACcomp[k] = c
     IUPACdeg[k] = d


# All IUPAC triplets that unambiguously translate to a single amino acid
codons = [
    'ATT','ATC','ATA',                   # Isoleucine (I)
    'CTT','CTC','CTA','CTG','TTA','TTG', # Leucine (L)
    'GTT','GTC','GTA','GTG',             # Valine (V)
    'TTT','TTC',                         # Phenylalanine (F)
    'ATG',                               # Methionine (M)
    'TGT','TGC',                         # Cysteine (C)
    'GCT','GCC','GCA','GCG',             # Alanine (A)
    'GGT','GGC','GGA','GGG',             # Glycine (G)
    'CCT','CCC','CCA','CCG',             # Proline (P)
    'ACT','ACC','ACA','ACG',             # Threonine (T)
    'TCT','TCC','TCA','TCG','AGT','AGC', # Serine (S)
    'TAT','TAC',                         # Tyrosine (Y)
    'TGG',                               # Tryptophan (W)
    'CAA','CAG',                         # Glutamine (Q)
    'AAT','AAC',                         # Asparagine (N)
    'CAT','CAC',                         # Histidine (H)
    'GAA','GAG',                         # Glutamic acid (E)
    'GAT','GAC',                         # Aspartic acid (D)
    'AAA','AAG',                         # Lysine (K)
    'CGT','CGC','CGA','CGG','AGA','AGG', # Arginine (R)
    'TAA','TAG','TGA'                    # Stop codon (-)
]

# single-letter IUPAC amino acids for each entry in codons
aminos = [
    'I','I','I',
    'L','L','L','L','L','L',
    'V','V','V','V',
    'F','F',
    'M',
    'C','C',
    'A','A','A','A',
    'G','G','G','G',
    'P','P','P','P',
    'T','T','T','T',
    'S','S','S','S','S','S',
    'Y','Y',
    'W',
    'Q','Q',
    'N','N',
    'H','H',
    'E','E',
    'D','D',
    'K','K',
    'R','R','R','R','R','R',
    '-','-','-'
]

# All codons with IUPAC ambiguities that code for a single amino acid
codons_ambiguous = [
    'ATY','ATW','ATM','ATD',             # Isoleucine (I)
    'CTM','CTR','CTW','CTS','CTY','CTK', # Leucine (L)
    'CTV','CTH','CTD','CTB','CTN',
    'YTA','YTG','TTR','YTR',
    'GTM','GTR','GTW','GTS','GTY','GTK', # Valine (V)
    'GTV','GTH','GTD','GTB','GTN',
    'TTY',                               # Phenylalanine (F)
    'TGY',                               # Cysteine (C)
    'GCM','GCR','GCW','GCS','GCY','GCK', # Alanine (A)
    'GCV','GCH','GCD','GCB','GCN',
    'GGM','GGR','GGW','GGS','GGY','GGK', # Glycine (G)
    'GGV','GGH','GGD','GGB','GGN',
    'CCM','CCR','CCW','CCS','CCY','CCK', # Proline (P)
    'CCV','CCH','CCD','CCB','CCN',
    'ACM','ACR','ACW','ACS','ACY','ACK', # Threonine (T)
    'ACV','ACH','ACD','ACB','ACN',
    'TCM','TCR','TCW','TCS','TCY','TCK', # Serine (S)
    'TCV','TCH','TCD','TCB','TCN','AGY',
    'TAY',                               # Tyrosine (Y)
    'CAR',                               # Glutamine (Q)
    'AAY',                               # Asparagine (N)
    'CAY',                               # Histidine (H)
    'GAR',                               # Glutamic acid (E)
    'GAY',                               # Aspartic acid (D)
    'AAR',                               # Lysine (K)
    'CGM','CGR','CGW','CGS','CGY','CGK', #Arginine (R)
    'CGV','CGH','CGD','CGB','CGN',
    'AGR','MGA','MGG','MGR',
    'TAR','TRA'                          # Stop codon (-)
]

# IUPAC amino acid that matches each element of codons_ambiguous
aminos_ambiguous = [
    'I','I','I','I',
    'L','L','L','L','L','L',
    'L','L','L','L','L',
    'L','L','L','L',
    'V','V','V','V','V','V',
    'V','V','V','V','V',
    'F',
    'C',
    'A','A','A','A','A','A',
    'A','A','A','A','A',
    'G','G','G','G','G','G',
    'G','G','G','G','G',
    'P','P','P','P','P','P',
    'P','P','P','P','P',
    'T','T','T','T','T','T',
    'T','T','T','T','T',
    'S','S','S','S','S','S',
    'S','S','S','S','S','S',
    'Y',
    'Q',
    'N',
    'H',
    'E',
    'D',
    'K',
    'R','R','R','R','R','R',
    'R','R','R','R','R',
    'R','R','R','R',
    '-','-'
]

quality_scores = [
    '!','"','#','$','%','&',
    "'",'(',')','*','+',
    ',','-','.','/','0',
    '1','2','3','4','5',
    '6','7','8','9',':',
    ';','<','=','>','?',
    '@','A','B','C','D',
    'E','F','G','H','I','J'
]

QUALhash = {}
for i in range(42):
    QUALhash[quality_scores[i]] = i

CODONhash = {}
for c,a in zip(codons, aminos):
    CODONhash[c] = a

AMBIGhash = {}
for c,a in zip(codons_ambiguous, aminos_ambiguous):
    AMBIGhash[c] = a

def which(x,value=True):
    """Returns a list of locations in x that satisfy value"""
    return [a for a,b in enumerate(x) if b==value]


def notwhich(x,value=0):
    """Returns a list of locations in x that do not satisty value"""
    return [a for a,b in enumerate(x) if b!=value]


def flatten(list_of_lists):
    """Collapses a list/tuple of lists into a single list"""
    return [item for sublist in list_of_lists for item in sublist]

def nuc_to_int(nuc_string):
    """Converts a string of IUPAC-encoded nucleotides
    to an integer array."""
    return [IUPACnum[i] for i in nuc_string]

def nuc_to_bin(nuc_string):
    """Converts a string of IUPAC-encoded nucleotides
    to an bytearray."""
    return bytearray([IUPACbin[ord(i)] for i in nuc_string])

def import_genome(genome_FASTA, split_on=' ', keep_case=True):
    """Reads FASTA file to a dict."""
    genome = {}
    chrom = 'none'
    genome_file = open(genome_FASTA)
    for line in genome_file:
        line = line.rstrip()
        if len(line)==0:
            continue
        if line[0] == '>':
            if chrom != 'none':
                genome[chrom] = ''.join(current_lines)
            chrom = line[1:].split(split_on)[0]
            current_lines = []
            continue
        
        if not keep_case:
            line = line.upper()
        current_lines.append(line)
    genome[chrom] = ''.join(current_lines)
    return genome

def number_chromosomes(genome_FASTA):
    """Returns a sorted index of chromosome starting positions in a genome."""
    genome = import_genome(genome_FASTA)
    
    running_count = 0
    chromosome_number = {}
    for c in sorted(genome.keys()):
        chromosome_number[c] = running_count
        running_count += 1
    
    return chromosome_number

def translate(codon):
    """Looks up a nucleotide triplet in the codon hashtables."""
    if len(codon) == 3:
        c = CODONhash.get(codon, '?')
    else:
        return ''
    
    if c == '?':
        return AMBIGhash.get(codon,'X')
    else:
        return c


def rc(sequence):
    """Returns reverse complement of a nucleotide string."""
    return ''.join(reversed([IUPACcomp.get(i,i) for i in sequence.upper()]))

def complement(sequence):
    """Returns complement of a nucleotide string."""
    return ''.join([IUPACcomp.get(i,i) for i in sequence.upper()])

def IUPAC_mismatches(sequence_a, sequence_b):
    """Calculates the number of mismatches present between
    two equal length strings a and b. IUPAC ambiguities do
    not count as mismatches, i.e. R matches A or G.
    """
    number_of_mismatches = 0
    a = sequence_a.upper()
    b = sequence_b.upper()
    for i,j in zip(a,b):
        if i == j:
            continue
        if i in IUPAChash[j] or j in IUPAChash[i]:
            continue
        number_of_mismatches += 1
    
    return number_of_mismatches


def longest_orf(sequence,allow_truncation=True):
    """Locates the longest open reading frame.
    
    Outputs a triple of:
        amino acid sequence (string)
        start:stop positions in nucleotide sequence (0-indexed int)
        Is the codon full or truncated? (bool)
    """
    sequence = sequence.upper()
    frame  = {}
    frame[0] = [sequence[i:(i+3)] for i in range(0,len(sequence),3)]
    frame[1] = [sequence[i:(i+3)] for i in range(1,len(sequence),3)]
    frame[2] = [sequence[i:(i+3)] for i in range(2,len(sequence),3)]
    orf = ''
    span = []
    stopless = False
    start_met = re.compile('^.*?(M.*)$')
    for f in [0,1,2]:
        translation = ''.join([translate(i) for i in frame[f]])
        potential_orfs = translation.split('-')
        stopless_list = [False]*(len(potential_orfs)-1)+[True]
        if not allow_truncation:
            potential_orfs = potential_orfs[:-1]
            stopless = stopless_list[:-1]
        for p,s in zip(potential_orfs,stopless_list):
            M_match = start_met.match(p)
            if M_match:
                o = M_match.groups()[0]
                if len(o) > len(orf):
                    orf = o
                    if s:
                        stopless = True
                        span = [i*3+f for i in re.search(
                            o,translation).span()]
                    else:
                        stopless = False
                        span = [i*3+f for i in re.search(
                            o+'-',translation).span()]
    return (orf,span,stopless)


def fasta_degeneracy(sequence):
    """Returns per-position level of ambiguity for a nucleotide string.
    
    Input:
        'AACAGRTCGNNNYGCHG'
    Output:
        [1,1,1,1,1,2,1,1,1,4,4,4,2,1,1,3,1]
    """
    return [IUPACdeg.get(i, 4) for i in sequence]
    

#def orf_mutations(orfsequence):
#    """Returns number of synonymous/nonsynonymous IUPAC ambiguities.
#    
#    Output a triple of ints:
#        # unambiguous nucleotides
#        # in-frame synonymous
#        # in-frame nonsynonymous
#    """
#    deg = fasta_degeneracy(orfsequence)
#    
#    all_codons = [orfsequence[i:(i+3)] for i in range(0,len(orfsequence),3)]
#    all_aa = [translate(i) for i in all_codons]


def to_regex(sequence):
    """Converts an IUPAC-formatted string to a regex string"""
    return ''.join([IUPAChash[i] for i in sequence.upper()])
