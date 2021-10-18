import numpy as np
import copy
from collections import namedtuple, Counter

bed_colors = {
    'U':'128,130,133',
    'C':'0,212,145',
    'S':'28,117,188',
    'E':'190,30,45',
    'SE':'109,74,116'
}

ELdata = namedtuple('ELdata', 'chrom source strand ranges splice s_tag e_tag capped weight')

class RNAseqMapping:
    def __init__(self, input_data):
        """Initializes a Read Object given a tuple of input data.
        Requires a chromosome, strand, source, weight, a sorted tuple of
        exon ranges and an array of booleans indicating which gaps between exons are splice junctions."""
        if input_data.source is None:
            self.source = None
        else:
            self.source = int(input_data.source)
        
        self.chrom = int(input_data.chrom)
        self.strand = input_data.strand
        self.ranges = sorted(input_data.ranges)
        self.splice = input_data.splice
        if self.strand == '.': # Terminal tag information is meaningless for nonstranded reads
            self.s_tag = self.e_tag = self.capped = False
        else:
            self.s_tag, self.e_tag, self.capped = input_data.s_tag, input_data.e_tag, input_data.capped
        
        self.weight = float(input_data.weight)
        self.span = (self.left(), self.right())
    
    def __eq__(self, other): return self.span == other.span
    def __gt__(self, other): return self.span > other.span
    def __ge__(self, other): return self.span >= other.span
    def __lt__(self, other): return self.span < other.span
    def __le__(self, other): return self.span <= other.span
    def __ne__(self, other): return self.span != other.span
    
    def __repr__(self):
        """Represents the read object with an ASCII character string:
            >>, <<, || represent plus, minus, and unstranded
            Ranges are connected by ^ (splice) or . (gap)"""
        if self.strand == '+':
            strandchar = '>'
        elif self.strand == '-':
            strandchar = '<'
        else:
            strandchar = '|'
        
        gapchar = ['^' if i else '_' for i in self.splice] + [strandchar]
        rangechar = ['{}-{}'.format(a,b) for a,b in self.ranges]
        return(''.join([a+b for a,b in zip(rangechar,gapchar)]))
    
    def left(self):
        return self.ranges[0][0]
    
    def right(self):
        return self.ranges[-1][-1]
    
    def get_length(self):
        """Returns the number of nucleotides covered by all blocks of the object."""
        return sum([b-a for a,b in self.ranges])
    
    def gaps(self):
        """Returns an array of 0-indexed (start, end) tuples of gaps between ranges"""
        if len(self.ranges) == 1:
            return []
        
        return [(self.ranges[i][-1], self.ranges[i+1][0]) for i in range(len(self.ranges)-1)]
    
    def junctions(self):
        """Returns an array of 0-indexed (start, end) tuples of intron locations"""
        j_array = []
        for i,j in enumerate(self.splice):
            if j:
                j_array += [(self.ranges[i][-1], self.ranges[i+1][0])]
        
        return j_array
    
    def overlaps(self, other):
        """Returns a boolean if the mapping range of self overlaps other."""
        l1, r1 = self.span
        l2, r2 = other.span
        if r1 < l2 or l1 > r2:
            return False
        
        return True
    
    def overlap_range(self, other):
        """Returns 0-indexed open range of overlap between self and other"""
        if not self.overlaps(other):
            return None
        
        edges = sorted([self.span, other.span])
        return (edges[1][0], edges[0][1])
    
    def ends_clash(self, other):
        """Returns a boolean of whether the combination of end tags between
        self and other can be substrings of the same object."""
        strand = sorted([self.strand, other.strand])
        if strand[0] == '+':
            assert strand[1] != '-'
            if self.s_tag and self.span[0] > other.span[0]: return True # Other is left of self's S tag
            if other.s_tag and other.span[0] > self.span[0]: return True # Self is left of other's S tag
            if self.e_tag and self.span[-1] < other.span[-1]: return True # Other is right of self's E tag
            if other.e_tag and other.span[-1] < self.span[-1]: return True # Self is right of other's E tag
        elif strand[0] == '-':
            assert strand[1] != '+'
            if self.s_tag and self.span[-1] < other.span[-1]: return True # Other is right of self's S tag
            if other.s_tag and other.span[-1] < self.span[-1]: return True # Self is right of other's S tag
            if self.e_tag and self.span[0] > other.span[0]: return True # Other is left of self's E tag
            if other.e_tag and other.span[0] > self.span[0]: return True # Self is left of other's E tag
        
        return False
    
    def is_compatible(self, other):
        """Self and other contain no attributes that demonstrate they could
        not be subsequences of a common longer molecule."""
        if self.chrom != other.chrom:
            return False
        
        if self.source != other.source:
            return False
        
        if self.strand != '.' and other.strand != '.' and self.strand != other.strand:
            return False
        
        # If self or other contain terminal tags, the tags must be compatible
        if self.ends_clash(other):
            return False
        
        overlap = self.overlap_range(other)
        if overlap is None:
            return True # No incompatibilities were found
        
        # If the two reads share a chrom and strand and overlap,
        # check the overlapping range for identical splice architecture
        j1 = [j for j in self.junctions() if j[0] > overlap[-1] and j[-1] < overlap[0]]
        j2 = [j for j in other.junctions() if j[0] > overlap[-1] and j[-1] < overlap[0]]
        if j1 == j2:
            return True
        
        return False
        
    def is_identical(self, other):
        """Boolean of whether two read objects share all nontrivial attributes"""
        if self.source != other.source: return False
        if self.chrom != other.chrom: return False
        if self.strand != other.strand: return False
        if self.ranges != other.ranges: return False
        if self.splice != other.splice: return False
        if self.s_tag != other.s_tag: return False
        if self.e_tag != other.e_tag: return False
        if self.capped != other.capped: return False
        return True
    
    def merge(self, other, check_compatibility=True):
        """Combines with another Mapping object. Must be compatible."""
        if check_compatibility:
            if not self.is_compatible(other):
                return False
        
        # Unify the strand information of the two objects
        self.s_tag = self.s_tag or other.s_tag
        self.e_tag = self.e_tag or other.e_tag
        self.capped = self.capped or other.capped
        self.weight = self.weight + other.weight
        if self.strand == '.':
            self.strand = other.strand
        
        o_range = self.overlap_range(other)
        if o_range is not None: # The two ranges overlap to some degree
            junctions = self.junctions() + other.junctions()
            self.ranges = collapse_blocks(self.ranges + other.ranges)
            new_gaps = self.gaps()
            self.splice = [gap in junctions for gap in new_gaps]
        elif self < other: # self is strictly left of other
            self.ranges = self.ranges + other.ranges
            self.splice = self.splice + [False] + other.splice # Join the two ranges with a gap
        else: # self is strictly right of other
            self.ranges = other.ranges + self.ranges
            self.splice = other.splice + [False] + self.splice # Join the two ranges with a gap
        
        self.span = (self.left(), self.right())
        return True
    
    def subtract(self, other):
        """Combines with another Mapping object. Must be compatible."""
        for r in other.ranges:
            self.ranges = subtract_range(self.ranges, r)
    
    def get_node_labels(self):
        """Returns a string with one label for each edge of each range in self.ranges."""
        if self.strand == '+':
            gapchar = 'DA'
            if self.s_tag:
                if self.capped:
                    startchar = 'C'
                else:
                    startchar = 'S'
            else:
                startchar = '.'
            
            if self.e_tag:
                endchar = 'E'
            else:
                endchar = '.'
        elif self.strand == '-':
            gapchar = 'AD'
            if self.s_tag:
                if self.capped:
                    endchar = 'C'
                else:
                    endchar = 'S'
            else:
                endchar = '.'
            
            if self.e_tag:
                startchar = 'E'
            else:
                startchar = '.'
        else:
            gapchar = '..'
            startchar = endchar = '.'
        
        return(''.join([startchar]+[gapchar if i else '..' for i in self.splice]+[endchar]))
    
    def write_as_elr(self, as_string=True):
        """Returns a string that represents the ReadObject
        in the end-labeled read (ELR) format"""
        block_ends = flatten(self.ranges)
        lengths = [block_ends[i]-block_ends[i-1] for i in range(1,len(block_ends))]
        labels = self.get_node_labels()
        EL_CIGAR = ''.join([str(a)+str(b) for a,b in zip(labels,lengths+[''])])
        elr_line = [self.chrom, self.left(), self.strand, EL_CIGAR, self.source, round(self.weight,2)]
        if as_string:
            return '\t'.join([str(i) for i in elr_line])
        else:
            return elr_line
     
    def write_as_bed(self, chrom_array, source_array, as_string=True, name='.'):
        """Returns a string that represents the ReadObject
        in a 15-column BED format"""
        labels = self.get_node_labels()
        l = labels[0]
        r = labels[-1]
        ends = l+r
        if ends in ['SE','ES']:
            rgb = bed_colors['SE']
        elif 'C' in ends:
            rgb = bed_colors['C']
        elif 'S' in ends:
            rgb = bed_colors['S']
        elif 'E' in ends:
            rgb = bed_colors['E']
        else:
            rgb = bed_colors['U']
        
        chromStart, blockStarts, blockSizes = explode_block_ranges(self.ranges)
        if self.source is None:
            source_str = ''
        else:
            source_str = source_array[self.source]
        
        bed_line = [
            chrom_array[self.chrom], chromStart, self.right(),
            name, '.', self.strand, 0, 0, rgb,
            len(self.ranges),
            ','.join([str(i-1) for i in blockSizes]),
            ','.join([str(i) for i in blockStarts]),
            self.weight, source_str, labels
        ]
        if as_string:
            return '\t'.join([str(i) for i in bed_line])
        else:
            return bed_line

class Path(RNAseqMapping):
    def __init__(self, input_data):
        RNAseqMapping.__init__(self, input_data)
    
    def get_read_depth(self):
        return self.weight/self.get_length()
    
    def is_complete(self, require_cap=False, maxgap=50):
        """Returns a boolean whether Path connects 
        a Start to an End with no gaps larger than maxgap"""
        if require_cap and not self.capped:
            return False
        
        if not self.s_tag or not self.e_tag:
            return False
        
        if any([(b-a)>maxgap for a,b in self.get_gaps()]):
            return False
        
        return True

class RNAseqDataset:
    def __init__(self, chrom_array=None, source_array=None, chrom_lengths=None, record_names=False):
        """Container for RNAseqMapping objects. Stores a reference dictionary for all
        chromosome names and sample names. Contains methods for parsing
        a variety of files into a collection of read objects."""
        self.read_list = []
        self.record_names = record_names
        if self.record_names:
            self.read_names = {}
        
        self.chrom_lengths = chrom_lengths
        
        self.chrom_dict = {}
        self.chrom_index = 0
        self.chrom_array = []
        if chrom_array is not None:
            for c in chrom_array:
                self.add_chrom(c)
        
        self.source_dict = {}
        self.source_index = 0
        self.source_array = []
        if source_array is not None:
            for s in source_array:
                self.add_source(s)
    
    def add_source(self, source_string):
        if source_string not in self.source_dict:
            self.source_array.append(source_string)
            self.source_dict[source_string] = self.source_index
            self.source_index += 1
    
    def add_chrom(self, chrom_string):
        if chrom_string not in self.chrom_dict:
            self.chrom_array.append(chrom_string)
            self.chrom_dict[chrom_string] = self.chrom_index
            self.chrom_index += 1
    
    def add_read_from_BED(self, bed_line, source_string=None, s_tag=False, e_tag=False, capped=False, gaps_are_junctions=False):
        input_data = parse_BED_line(bed_line, self.chrom_dict, self.source_dict, source_string, s_tag, e_tag, capped, gaps_are_junctions)
        
        if type(input_data.chrom) is str: # Chromosome wasn't in chrom_dict
            chrom_string = input_data.chrom
            input_data = input_data._replace(chrom=self.chrom_index)
            self.add_chrom(chrom_string)
        
        if type(input_data.source) is str: # Source wasn't in source_dict
            source_string = input_data.source
            input_data = input_data._replace(source=self.source_index)
            self.add_source(source_string)
        
        new_read = RNAseqMapping(input_data)
        self.read_list.append(new_read)
        if self.record_names:
            name = bed_line.split('\t')[3].split('.')[0].upper()
            if name in self.read_names:
                self.read_names[name].merge(self.read_list[-1],check_compatibility=False)
            else:
                self.read_names[name] = self.read_list[-1]
    
    def add_read(self, input_data):
        if type(input_data.chrom) is str: # Chromosome wasn't in chrom_dict
            chrom_string = input_data.chrom
            input_data = input_data._replace(chrom=self.chrom_index)
            self.add_chrom(chrom_string)
        
        if type(input_data.source) is str: # Source wasn't in source_dict
            source_string = input_data.source
            input_data = input_data._replace(source=self.source_index)
            self.add_source(source_string)
        
        new_read = RNAseqMapping(input_data)
        self.read_list.append(new_read)
    
    def add_read_from_ELR(self, elr_line):
        new_read = elr_to_readobject(elr_line)
        self.read_list.append(new_read)
    
    def pop_read(self, read_format='elr', as_string=True):
        """Remove the last read added to the stack and write it in 'format'.
        """
        if read_format.lower() == 'elr':
            return(self.read_list.pop().write_as_elr(as_string))
        elif read_format.lower() == 'bed':
            return(self.read_list.pop().write_as_bed(self.chrom_array, self.source_array, as_string))
    
    def dump_header(self):
        """Returns an array of strings that describe chrom_dict and source_dict of the Dataset."""
        header_list = []
        for i,c in enumerate(self.chrom_array):
            header_list += ['#C {} {}'.format(i, c)]
        
        for i,s in enumerate(self.source_array):
            header_list += ['#S {} {}'.format(i, s)]
        
        return header_list




bp_typeorder = {'SOURCE':-1, 'E':0, 'A':1, 'D':2, 'S':3, 'SINK':4} # Sort order for branchpoint types
class BranchPoint:
    """Represents a reference point for a Start, End, Donor, or Acceptor site."""
    def __init__(self, branchtype, strand, pos, weight):
        self.strand = strand
        assert self.strand in ['+','-']
        self.branchtype = branchtype
        self.pos = pos
        self.weight = float(weight)
        if self.strand == '+':
            self.comparator = (self.pos, bp_typeorder[self.branchtype])
        else:
            self.comparator = (self.pos, -bp_typeorder[self.branchtype])
    
    def __repr__(self):
        return '{}{}{} ({})'.format(self.branchtype, self.strand, self.pos, self.weight)
    
    def __eq__(self, other): return self.comparator == other.comparator
    def __ne__(self, other): return self.comparator != other.comparator
    def __gt__(self, other): return self.comparator >  other.comparator
    def __ge__(self, other): return self.comparator >= other.comparator
    def __lt__(self, other): return self.comparator <  other.comparator
    def __le__(self, other): return self.comparator <= other.comparator

class StartPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight, extend=50):
        BranchPoint.__init__(self, 'S', strand, pos, weight)
        self.extend = extend
        self.left = self.right = self.pos
        self.range = range(self.left-self.extend, self.right+self.extend+1)
    
    def merge(self, other):
        """Merge other into self"""
        self.weight += other.weight
        self.left = min(self.left, other.left)
        self.right = max(self.right, other.right)
        self.range = range(self.left-self.extend, self.right+self.extend+1)

class EndPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight, extend=50):
        BranchPoint.__init__(self, 'E', strand, pos, weight)
        self.extend = extend
        self.left = self.right = self.pos
        self.range = range(self.left-self.extend, self.right+self.extend+1)
    
    def merge(self, other):
        """Merge other into self"""
        self.weight += other.weight
        self.left = min(self.left, other.left)
        self.right = max(self.right, other.right)
        self.range = range(self.left-self.extend, self.right+self.extend+1)

class DonorPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight):
        BranchPoint.__init__(self, 'D', strand, pos, weight)

class AcceptorPoint(BranchPoint):
    """Represents a transcript start site"""
    def __init__(self, strand, pos, weight):
        BranchPoint.__init__(self, 'A', strand, pos, weight)

class SourcePoint(BranchPoint):
    """Special class of branchpoint, represents global source.
    Can only form edges with S edges"""
    def __init__(self, strand, weight):
        BranchPoint.__init__(self, 'SOURCE', strand, -float('inf'), weight)

class SinkPoint(BranchPoint):
    """Special class of branchpoint, represents global source.
    Can only form edges with S edges"""
    def __init__(self, strand, weight):
        BranchPoint.__init__(self, 'SINK', strand, float('inf'), weight)

Capacity = namedtuple('Capacity', 'plus minus ns')
class Vertex:
    """Base class for a vertex in a locus graph object"""
    def __init__(self):
        self.capacity = Capacity(0, 0, 0)
    
    def __eq__(self, other): return self.comparator == other.comparator
    def __ne__(self, other): return self.comparator != other.comparator
    def __gt__(self, other): return self.comparator >  other.comparator
    def __ge__(self, other): return self.comparator >= other.comparator
    def __lt__(self, other): return self.comparator <  other.comparator
    def __le__(self, other): return self.comparator <= other.comparator

class TerminalVertex(Vertex):
    """Child class for the external vertices of a locus."""
    def __init__(self, is_start=False):
        Vertex.__init__(self)
        self.capped = float(0)
        if is_start:
            self.side = 'L'
            self.comparator = -float('inf')
        else:
            self.side = 'R'
            self.comparator = float('inf')
    
    def __repr__(self):
        return '<TV {}>'.format(self.side)

class InternalVertex(Vertex):
    """Child class for an internal vertex of a locus graph object"""
    def __init__(self, left, right):
        Vertex.__init__(self)
        self.left = left
        self.right = right
        self.range = range(self.left, self.right)
        self.gaps = np.zeros((self.right-self.left), dtype=np.bool) # Records whether each position is covered by >=1 read
        self.comparator = self.left
    
    def __repr__(self):
        return '<IV {}-{}>'.format(self.left,self.right)

class Edge:
    """Base class for an edge in a locus graph object"""
    def __init__(self):
        pass

class Adjacency(Edge):
    """An edge connecting two adjacent vertices, not associated
    with any branchpoints. Carries +, -, and . capacity"""
    def __init__(self, left_vertex, right_vertex):
        Vertex.__init__(self)

class Bridge(Edge):
    """An edge connecting two vertices via one (S, E), or two (D|A, A|D) branchpoints.
    Carries capacity that matches the strand of the branchpoint(s)"""
    def __init__(self, left_, right):
        Vertex.__init__(self)

class OrderedArray:
    """Maintains an array in sort order and a lookup dict for the location of each
    item that was added."""
    def __init__(self, array=[]):
        self.n = len(array)
        self.ordered_array = sorted([(v,i) for i,v in enumerate(array)])
    
    def __repr__(self):
        return 'OA({})'.format(self.n)
    
    def ordered_index(self, position=None):
        if position is None:
            return [i for v,i in self.ordered_array]
        else:
            return self.ordered_array[position][-1]
    
    def ordered_values(self, position = None):
        if position is None:
            return [v for v,i in self.ordered_array]
        else:
            return self.ordered_array[position][0]
    
    def __iter__(self):
        return iter(self.ordered_values())
    
    def __len__(self):
        return len(self.ordered_array)
    
    def __add__(self,other):
        new_OA = copy.deepcopy(self)
        for item in other:
            new_OA.add(item)
        
        return new_OA
    
    def get_insert_pos(self, search_value, init=None):
        """Runs a binary search, returns an insert position where 'search_value'
        fits in array. Assumes array to be sorted. In the case of a tie,
        priority goes to the existing item."""
        insert_pos = -1
        if init is None:
            first = 0
            last = len(self.ordered_array) - 1
        else: # Start at an initialized index position (if a good guess of the correct index exists)
            if search_value == self.ordered_array[init]: # Exact match, insert pos is right of init
                return init + 1
            elif search_value > self.ordered_array[init]: # insert pos must be right of init
                first = init + 1
                last = len(self.ordered_array) - 1
            else: # insert pos must be left of init
                first = 0
                last = init - 1
        
        while (first <= last) and insert_pos == -1:
            mid = (first + last)//2 # Get the midpoint of the search space
            if search_value == self.ordered_array[mid]: # Exact match, insert pos is right of mid
                insert_pos = mid + 1
            elif search_value > self.ordered_array[mid]: # insert pos must be right of the midpoint
                first = mid + 1
            else: # insert pos must be left of the midpoint
                last = mid - 1
        
        if insert_pos == -1: # An exact match wasn't found, but the search converged
            insert_pos = first
        
        return insert_pos
    
    def add(self, item, init=None):
        """Add an item to the ordered array."""
        self.n += 1
        ipos = self.get_insert_pos((item, self.n), init)
        self.ordered_array.insert(ipos, (item, self.n))
    
    def delete(self, index):
        """Deletes the item at index position in the ordered array."""
        self.n -= 1
        del self.ordered_array[index]
    
    def add_list(self, list_to_add):
        for item in list_to_add:
            self.add(item)
    
    def probe(self, item, init=None):
        """Return the insert position of an item if it were to be added to the array"""
        return self.get_insert_pos((item, self.n + 1), init)


class LocusGraph:
    def __init__(self, list_of_reads=[]):
        self.leftmost = self.rightmost = None
        self.reads = tuple(list_of_reads) # Cannot be mutated
        self.weight = float(0)
        
        self.bp = {}
        if len(self.reads) > 0:
            self.generate_branchpoints()
            self.build_graph()
        else:
            self.paths = OrderedArray()
            self.bp['+'] = OrderedArray() # Empty array for storing an ordered list of plus-stranded branchpoints
            self.bp['-'] = OrderedArray() # Empty array for storing an ordered list of minus-stranded branchpoints
    
    def distance_from_edge(self, pos, strand, type):
        """ Given a position (int), strand (+,-), and type (S,E,D,A),
        returns the the distance to the most extreme edge of the locus."""
        if type in ['S','D']: # Consider upstream
            if strand == '+':
                return pos - self.leftmost
            else:
                return self.rightmost - pos
        else: # Consider downstream
            if strand == '+':
                return self.rightmost - pos
            else:
                return pos - self.leftmost
    
    def generate_branchpoints(self):
        """Makes a collection of positions in the locus that act as dividing
        points for all nodes of the graph. Every donor (D) and acceptor (A) 
        site are retained as-is, but start (S) and end (E) are collapsed into
        a set of reference points."""
        bp_dict = dict()
        for branchtype in ['S','E','D','A']:
            for strand in ['+','-']:
                bp_dict[branchtype+strand] = Counter()
        
        bp_dict['SOURCE'] = Counter()
        bp_dict['SINK'] = Counter()
        
        self.leftmost = self.reads[0].span[0]
        self.rightmost = self.reads[-1].span[-1]
        for read in self.reads:
            if read.span[-1] > self.rightmost:
                self.rightmost = read.span[-1]
            
            self.weight += read.weight
            if read.strand == '+':
                for l,r in read.junctions():
                    bp_dict['D+'][l] += read.weight
                    bp_dict['A+'][r] += read.weight
                
                if read.s_tag:
                    pos = read.span[0]
                    bp_dict['S+'][pos] += read.weight
                    bp_dict['SOURCE']['+'] += read.weight
                
                if read.e_tag:
                    pos = read.span[-1]
                    bp_dict['E+'][pos] += read.weight
                    bp_dict['SINK']['+'] += read.weight
            elif read.strand == '-':
                for l,r in read.junctions():
                    bp_dict['A-'][l] += read.weight
                    bp_dict['D-'][r] += read.weight
                
                if read.e_tag:
                    pos = read.span[0]
                    bp_dict['E-'][pos] += read.weight
                    bp_dict['SINK']['-'] += read.weight
                
                if read.s_tag:
                    pos = read.span[-1]
                    bp_dict['S-'][pos] += read.weight
                    bp_dict['SOURCE']['-'] += read.weight
        
        # Populate the two BP arrays with D/A sites
        self.bp['+'] = OrderedArray(
            [DonorPoint('+', k, v) for k,v in bp_dict['D+'].items()] + 
            [AcceptorPoint('+', k, v) for k,v in bp_dict['A+'].items()]
        )
        self.bp['-'] = OrderedArray(
            [DonorPoint('-', k, v) for k,v in bp_dict['D-'].items()] + 
            [AcceptorPoint('-', k, v) for k,v in bp_dict['A-'].items()]
        )
        
        # Generate a priority queue for adding S/E branchpoints
        for strand in ['+','-']:
            start_order = [(i, self.distance_from_edge(v[0], strand, 'S'), -v[1], 'S', v[0]) for i,v in enumerate(bp_dict['S'+strand].most_common())]
            end_order = [(i, self.distance_from_edge(v[0], strand, 'E'), -v[1], 'E', v[0]) for i,v in enumerate(bp_dict['E'+strand].most_common())]
            priority_ranked_SE_array = sorted(start_order + end_order)
            for o,d,w,t,p in priority_ranked_SE_array:
                weight = -w
                if t == 'S':
                    item = StartPoint(strand, p, weight)
                else:
                    item = EndPoint(strand, p, weight)
                
                self.add_endpoint(item)
        
        terminal_branchpoints = OrderedArray([
            SourcePoint('+',bp_dict['SOURCE']['+']),SourcePoint('-',bp_dict['SOURCE']['-']), 
            SinkPoint('+',bp_dict['SINK']['+']), SinkPoint('-',bp_dict['SINK']['-'])
        ])
        # Collapse all branchpoints into a single fixed tuple
        self.branchpoints = tuple(self.bp['+'] + self.bp['-'] + terminal_branchpoints)
        self.D = {}
        self.D['+'] = dict([(b.pos,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'D' and b.strand == '+'])
        self.D['-'] = dict([(b.pos,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'D' and b.strand == '-'])
        self.A = {}
        self.A['+'] = dict([(b.pos,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'A' and b.strand == '+'])
        self.A['-'] = dict([(b.pos,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'A' and b.strand == '-'])
        
        # Make a constant-time lookup to point to the appropriate branchpoint for any 
        # Read through all S/E positions and make a constant-time lookup table for which branchpoint each belongs to
        S_plus = [(b,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'S' and b.strand == '+']
        S_minus = [(b,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'S' and b.strand == '-']
        E_plus = [(b,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'E' and b.strand == '+']
        E_minus = [(b,i) for i,b in enumerate(self.branchpoints) if b.branchtype == 'E' and b.strand == '-']
        self.S = {'+':{},'-':{}}
        self.E = {'+':{},'-':{}}
        for k in bp_dict['S+'].keys():
            for b,i in S_plus:
                if k in b.range:
                    self.S['+'][k] = i
                    break
        
        for k in bp_dict['S-'].keys():
            for b,i in S_minus:
                if k in b.range:
                    self.S['-'][k] = i
                    break
        
        for k in bp_dict['E+'].keys():
            for b,i in E_plus:
                if k in b.range:
                    self.E['+'][k] = i
                    break
        
        for k in bp_dict['E-'].keys():
            for b,i in E_minus:
                if k in b.range:
                    self.E['-'][k] = i
                    break
    
    def add_endpoint(self, item):
        """Add an S/E BranchPoint to the bp dict of OrderedArrays.
        If at least one same-type branchpoint is in range, merge item
        with the higher-priority branchpoint. If item is flanked by two
        neighboring branchpoints, merge both item and secondary into the first."""
        array = self.bp[item.strand]
        insert_site = array.probe(item)
        
        left, lp = array.ordered_array[insert_site-1] # Get the left neighbor and its priority order
        if insert_site < array.n:
            right, rp = array.ordered_array[insert_site] # Get the right neighbor and its priority order
        else: # No right-sided neighbor exists
            right = rp = None
        
        l_in_range = r_in_range = False
        if left.branchtype == item.branchtype:
            if item.pos in left.range: # Item could be merged with left
                l_in_range = True
        
        if right is not None:
            if right.branchtype == item.branchtype:
                if item.pos in right.range: # Item could be merged with right
                    r_in_range = True
        
        if l_in_range:
            if r_in_range:
                if lp < rp: # left is prioritized over right
                    left.merge(item)
                    left.merge(right)
                    array.delete(insert_site+1) # Remove right from the array
                else: # right is prioritized over left
                    right.merge(item)
                    right.merge(left)
                    array.delete(insert_site) # Remove left from the array
            else: # only left can merge item
                left.merge(item)
        elif r_in_range: # only right can merge item
            right.merge(item)
        else: # no compatible neighbors
            array.add(item)
    
    def remove_branchpoint(self, strand, index):
        """Remove the BranchPoint at index in the strand bp dict.
        If the neighboring branchpoints are the same type and within range of each other,
        merge the smaller into the larger."""
        array = self.bp[strand]
        array.delete(index)
        
        left, lp = array.ordered_array[index-1] # Get the left neighbor and its priority order
        if index < array.n:
            right, rp = array.ordered_array[index] # Get the right neighbor and its priority order
            if left.branchtype == right.branchtype: # Both neighbors were the same type
                if left.branchtype in ['S','E']: # Branchtype is terminal
                    if left.range[-1] >= right.range[0]: # Neighbor ranges overlap
                        if lp < rp: # left is prioritized over right
                            left.merge(right)
                            array.delete(index) # Remove right from the array
                        else: # right is prioritized over left
                            right.merge(left)
                            array.delete(index-1) # Remove left from the array
    
    def update_edge(self, source_index, sink_index, weight):
        """Adds a new edge or increments an existing edge
        given two branchpoints.
        """
        source = self.branchpoints[source_index]
        sink = self.branchpoints[sink_index]
        assert source.strand == sink.strand
        
    
    def build_graph(self, min_overhang=1):
        """After branchpoints are identified, populate two tables that store information about
        whether every branchpoint and every vertex, respectively, is (1) present in the read,
        (-1) incompatible with the read, or (0) not measured by the read.
        Construct all Vertices and Edges in one pass through the reads."""
        bp_positions = sorted(list(set([i.pos for i in self.branchpoints])))[1:-1]
        self.blocks = tuple((a,b) for a,b in zip(bp_positions[:-1], bp_positions[1:])) # Immutable after init, makes a vertex between all pairs of nonterminal branchpoint positions
        self.branchpoint_table = np.zeros((len(self.reads), len(self.branchpoints)), dtype=np.int8) # Container for membership of each branchpoint in each read (-1 False, 1 True, 0 Unmeasured)
        self.branchpoint_strand = np.array([1 if i.strand=='+' else 0 for i in self.branchpoints], dtype=np.bool) # Boolean array for fast lookup of branchpoint strands
        # The first two and last two indices are ALWAYS: [SOURCE+, SOURCE-,..., SINK-, SINK+]
        self.block_table = np.zeros((len(self.reads), len(self.blocks)), dtype=np.int8) # Container for membership of each vertex in each read (-1 False, 1 True, 0 Unmeasured)
        self.vertices = [InternalVertex(l,r) for l,r in self.blocks]+[TerminalVertex(is_start=True), TerminalVertex(is_start=False)]
        self.graph = Graph(len(self.vertices)) # Vertices at position -2 and -1 are reserved for the left and right Terminal Vertex, respectively
        
        left_index = 0
        right_index = 0
        li_start = self.blocks[left_index][0]
        ri_start = self.blocks[right_index][0]
        for i,read in enumerate(self.reads): # Read through the reads once, cataloging blocks and branchpoints present/absent
            if left_index >= len(self.blocks):
                break
            
            if read.strand == '+': # Plus-stranded read
                self.branchpoint_table[i, self.branchpoint_strand == False] = -1 # No minus-stranded branchpoints can be part of this read
                if read.s_tag:
                    source_index,sink_index = 0, self.S['+'][read.span[0]]
                    self.branchpoint_table[i,[source_index,sink_index]] = 1 # Update membership of SOURCE+ and the positional branchpoint
                    self.update_edge(source_index, sink_index, read.weight)
                if read.e_tag:
                    source_index, sink_index = self.E['+'][read.span[-1]], -1
                    self.branchpoint_table[i,[source_index, sink_index]] = 1 # Update membership of SINK+ and the positional branchpoint
                    self.update_edge(source_index, sink_index, read.weight)
                
                for d,a in read.junctions():
                    source_index, sink_index = self.D['+'][d], self.A['+'][a]
                    self.branchpoint_table[i, [source_index, sink_index]] = 1
                    self.update_edge(source_index, sink_index, read.weight)
                    
            elif read.strand == '-': # Minus-stranded start read
                self.branchpoint_table[i, self.branchpoint_strand == True] = -1 # No plus-stranded branchpoints can be part of this read
                if read.s_tag:
                    source_index, sink_index = 1, self.S['-'][read.span[-1]]
                    self.branchpoint_table[i,[source_index, sink_index]] = 1 # Update membership of SOURCE- and the positional branchpoint
                    self.update_edge(source_index, sink_index, read.weight)
                
                if read.e_tag:
                    source_index, sink_index = self.E['-'][read.span[0]], -2
                    self.branchpoint_table[i,[source_index, sink_index]] = 1 # Update membership of SINK- and the positional branchpoint
                    self.update_edge(source_index, sink_index, read.weight)
                
                for a,d in read.junctions():
                    source_index, sink_index = self.D['-'][d], self.A['-'][a]
                    self.branchpoint_table[i,[source_index, sink_index]] = 1
                    self.update_edge(source_index, sink_index, read.weight)
            
            # for i,r in enumerate(read.ranges): # Run through each block range of read
                # li = ri = left_index
                # while l > self.blocks[li][-1]:
                    # li_start = self.blocks[li][0]
                    # li += 1
                    # if li >= len(self.blocks):
                        # break
                
                # while r > self.blocks[ri][-1]:
                    # ri_start = self.blocks[ri][0]
                    # ri += 1
                    # if ri >= len(self.blocks):
                        # break
                
                # if li >= len(self.blocks):
                        # break
                
                # print(l, r, li_start, ri_start, li, ri)
                # self.block_gaps[li][(l-li_start):(r-li_start)] = 1
                # if ri != li:
                    # self.block_gaps[ri][(l-ri_start):(r-ri_start)] = 1
    
        
def longest_gap(np_array):
    """Given a numpy array, returns the length of the longest chain of 'falsy'
    values in the array"""
    transitions = list(np.where(np.diff(np_array==0))[0])
    if not np_array[0]:
        transitions.insert(0, -1)
    
    if not np_array[-1]:
        transitions.append(len(np_array)-1)
    
    return max([b-a for a,b in zip(transitions[:-1:2],transitions[1::2])])

# ENVIRONMENT SETUP: CLASSES
# Rewritten in Python from a Java implementation in 'Algorithms, 4th edition' by Robert Sedgewick, Kevin Wayne
class Graph():
    """An implementation of the Graph() class that allows
    parallel edges by storing edges as a list."""
    def __init__(self, V=0):
        self.directed = True
        self.V = V
        self.E = 0
        self.adj = {}
        self.vertices = list(range(self.V))
        for v in self.vertices:
            self.adj[v] = []
        
        self.edges = []
    
    def __repr__(self):
        return '<V={}, E={}>'.format(
            self.V,
            self.E
        )
    
    def add_edge(self, v, w, directed, weight=None):
        e = fEdge(v, w, weight)
        self.edges.append(e)
        self.adj[v].append(e)
        self.adj[w].append(e) # fEdge added to the reverse adjacency list so backward flow can be used
        self.E += 1
    
    def show_edges(self):
        if self.directed:
            connector = '->'
        else:
            connector = '--'
        
        for v in self.adj.keys():
            for e in self.adj[v]:
                # Get the other vertex of each Edge object in adj[v]
                w = e.other(v)
                print('{}{}{}'.format(v,connector,w))
    
    def summary(self):
        print('Max Degree: {}\nAverage Degree: {}\nSelf Loops: {}'.format(
            self.maxDegree(),
            self.averageDegree(),
            self.numberOfSelfLoops()
        )) 
    
    def degree(self, v):
        """Returns the degree of vertex v."""
        return len(self.adj[v])

    def maxDegree(self):
        """Returns the highest degree of any vertex in the graph."""
        max_dv = 0
        for v in self.vertices:
            dv = self.degree(v)
            if dv > max_dv:
                max_dv = dv
        
        return max_dv

    def averageDegree(self):
        """Returns the mean connectivity of the graph."""
        if self.directed:
            return float(self.E / self.V)
        else:
            return float(2.0 * self.E / self.V)

    def numberOfSelfLoops(self):
        """Returns the number of instances where 
        both ends of an edge are the same vertex."""
        loop_count = 0
        for v in self.vertices:
            # Get other vertext of Edge object and see if it is the same
            loop_count += sum([1 for e in self.adj[v] if e.other(v) == v])
        
        if self.directed:
            return loop_count
        else:
            return loop_count/2
    
    def reverse(self):
        """(For directed graphs) Returns a graph where all
        edges are reversed."""
        if not self.directed:
            return self
        
        rGraph = Graph(directed=self.directed)
        rGraph.V = self.V
        rGraph.E = self.E
        rGraph.vertices = set(range(rGraph.V))
        rGraph.adj = {}
        for v in rGraph.vertices:
            rGraph.adj[v] = []
        
        for v in self.adj.keys():
            for e in self.adj[v]:
                # Get edge properties
                w = e.other(v)
                weight = e.weight
                rGraph.add_edge(w, v, self.directed, weight)
        return rGraph


# Flow edges are directed and carry both a capacity and a flow
class fEdge():
    """A class for directed edges to be represented in a Graph()."""
    def __init__(self, v, w, capacity=0):
        self.v = int(v)
        self.w = int(w)
        self.flow = 0
        self.capacity = capacity
        self.directed = True
    
    def __eq__(self, e): return self.capacity == e.capacity
    def __ne__(self, e): return self.capacity != e.capacity
    def __gt__(self, e): return self.capacity > e.capacity
    def __lt__(self, e): return self.capacity < e.capacity
    def __le__(self, e): return self.capacity <= e.capacity
    def __ge__(self, e): return self.capacity >= e.capacity
    def __add__(self, e): return self.capacity + e.capacity
    def __sub__(self, e): return self.capacity - e.capacity
    def __neg__(self): return fEdge(self.v, self.w, -self.capacity, self.directed)
    
    def __repr__(self):
        separator = '->'
        return '<{}{}{} ({}/{})>'.format(
            self.v,
            separator,
            self.w,
            self.flow,
            self.capacity
        )
    
    def residualCapacityTo(self, v):
        """Returns the residual forward or reverse flow
        available to vertex v"""
        if v == self.w:
            # Flow is forward
            return self.capacity - self.flow
        else:
            # Reverse flow
            return self.flow
    
    def addResidualFlowTo(self, v, delta):
        """Changes self.flow in direction v by amount delta."""
        if v == self.w:
            # Flow is forward
            self.flow += delta
        else:
            # Reverse flow
            self.flow += -delta
        
        # Sanity check: Flow can never exceed capacity
        assert self.flow <= self.capacity and self.flow >= 0
        
    def source(self):
        """Returns the source of the edge"""
        return self.v
    
    def sink(self):
        """Returns the sink of the edge"""
        return self.w
    
    def other(self, v):
        """Returns the other vertex of the edge"""
        if v == self.v:
            return self.w
        else:
            return self.v
    
    def compareTo(self, other):
        """Compares the capacity of this Edge object to another Edge object.
        Returns 'less than' (-1), 'greater than' (+1), or 'equal' (0)."""
        if self.capacity < other.capacity:
            return -1
        elif self.capacity > other.capacity:
            return 1
        else:
            return 0

class ffMaxflow():
    """Calculates Maxflow from source s to target t using
    the Ford-Fulkerson algorithm."""
    def __init__(self, G, s, t):
        self.G = copy.deepcopy(G)
        self.vertices = self.G.vertices
        self.s = s
        self.t = t
        self.value = 0 # Value of the flow
        self.marked = [False]*self.G.V
        self.edgeTo = [None]*self.G.V
        self.FordFulkerson(self.s,self.t)
    
    def hasAugmentingPath(self):
        """Determines using breadth-first search whether the network has a path
        that would increase the total flow value from s to t"""
        self.marked = [False]*self.G.V
        self.edgeTo = [None]*self.G.V
        q = [] # Initialize a queue for BFS
        q.append(self.s) # Add source vertex to the queue
        self.marked[self.s] = True
        while len(q) > 0:
            v = q.pop(0)
            for e in self.G.adj[v]: # For every edge in the adjacency table
                w = e.other(v)
                if e.residualCapacityTo(w) > 0 and not self.marked[w]:
                    # Capacity exists in the edge and the sink vertex is not yet part of the cut
                    self.edgeTo[w] = e # Connect w to the cut by e
                    self.marked[w] = True
                    q.append(w)
        
        return self.marked[self.t] # Can t be reached by s via the network of residuals?
    
    def inCut(self, v):
        return self.marked[v]
    
    def FordFulkerson(self, s, t):
        while self.hasAugmentingPath(): # As long as a path that increases flow exists...
            # Determine the bottleneck capacity
            bottle = float('inf')
            v = t
            while v != s:
                bottle = min([bottle, self.edgeTo[v].residualCapacityTo(v)]) # Reduce the bottleneck value
                v = self.edgeTo[v].other(v) # Run backwards up the path to s
            
            # Augment the flow
            v = t
            while v != s:
                self.edgeTo[v].addResidualFlowTo(v, bottle)
                v = self.edgeTo[v].other(v) # Run backwards up the path to s
            
            self.value += bottle # Add the bottleneck value to the total flow f

class RNAseqMultimapper:
    """Contains multiple RNAseqMapping instances.
    The mapping_list is sorted so that an identical read will have an idential mapping_list."""
    def __init__(self, mapping_list):
        self.mapping_list = sorted(mapping_list) # Sort array of mappings by (start, end) location
        self.span = self.mapping_list[0].span # Index by lowest mapping coordinates
    
    def __eq__(self, other): return self.span == other.span
    def __gt__(self, other): return self.span > other.span
    def __ge__(self, other): return self.span >= other.span
    def __lt__(self, other): return self.span < other.span
    def __le__(self, other): return self.span <= other.span
    def __ne__(self, other): return self.span != other.span

class ReadStack():
    '''Temporary stack of reads for compiling identical species'''
    def __init__(self,chrom=None,pos=None,strand=None,species={}):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand
        self.species = species
    
    def __repr__(self):
        return '<{}:{}>'.format(
            self.chrom,
            self.pos
        )
    
    def populate(self,strand,endpos,string,value):
        '''Add a read to the species dict'''
        if strand not in self.species:
            self.species[strand] = {}
        if endpos not in self.species[strand]:
            self.species[strand][endpos] = {}
        self.species[strand][endpos][string] = \
            self.species[strand][endpos].get(string, float(0)) + value
    
    def dump(self,counter,digits = 3):
        '''Dump the species dict to stdout in a sorted BED format'''
        new_counter = counter.copy()
        # Print all species as BED lines
        for strand in self.species.keys():
            entries = sorted([(k,v) for k,v in self.species[strand].items()])
            for endpos,posdict in entries:
                species = sorted([(k,v) for k,v in posdict.items()])
                if args.BED12:
                    # Output an EndMap BED12 formatted line for each species
                    for string,score in species:
                        label,rgb,blocknum,blocksizes,blockstarts,samplename = string.split('\t')
                        new_counter[strand] += 1
                        out_string = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            self.chrom,self.pos,endpos,
                            label,1000,strand,
                            0,0,rgb,blocknum,blocksizes,blockstarts,
                            round(score,digits),samplename,label
                        )
                        print(out_string)
                else:
                    for string,score in species:
                        head,sequence,tail,length = string.split('\t')
                        new_counter[strand] += 1
                        print(
                            '\t'.join(
                                [
                                    str(i) for i in [
                                        self.chrom,
                                        self.pos,
                                        endpos,
                                        'p.{}'.format(new_plus),
                                        round(score,digits),
                                        strand,
                                        head,
                                        sequence,
                                        tail,
                                        length
                                    ]
                                ]
                            )
                        )
        
        return new_counter

#######################################################
def array_to_blocks(array):
    """Breaks an array into a collections of (start,end) INCLUSIVE doubles
    that define the coordinates of all contiguous blocks in the array"""
    clean_array = sorted(list(set([i for i in array if type(i) is int])))
    if len(clean_array) == 0:
        raise StopIteration
    
    start = clean_array[0]
    end = clean_array[0]
    if len(array) == 1:
        yield (start,end)
        raise StopIteration
    
    for i in clean_array[1:]:
        if i-1 == end:
            end = i
        else:
            yield (start,end)
            start = i
            end = i
    
    yield (start,end)

def overlap_type(range_a,range_b):
    """Returns a list of relationships
    between the edges of range_a and range_b"""
    return [
        int(range_a[0] >= range_b[0]),
        int(range_a[0] >= range_b[1]),
        int(range_a[1] >= range_b[0]),
        int(range_a[1] >= range_b[1])
    ]

def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

def collapse_blocks(list_of_doubles):
    """Iterates over a sorted list of 0-indexed open (start,end) doubles
    and merges those that overlap to return a set of nonoverlapping blocks"""
    list_of_doubles.sort()
    new_blocks = []
    last_block = None
    for current_block in list_of_doubles:
        if last_block:
            # Get the overlap type of the two adjacent blocks
            overlap = sum(overlap_type(current_block,last_block))
            if overlap > 0 and overlap < 4:
                # The two blocks do overlap to some degree
                current_block = (min(last_block[0],current_block[0]),max(last_block[-1],current_block[-1]))
            else:
                new_blocks += [last_block]
        
        # Update the last block
        last_block = current_block
    
    # Resolve the last block
    new_blocks += [last_block]
    return new_blocks

def get_block_ranges(chromStart, blockStarts, blockSizes):
    """Reorganizes the BED12 block format into a list of 0-indexed open (start,end) doubles."""
    if type(blockSizes) is str:
        sizes = [int(i) for i in blockSizes.split(',')]
        blockSizes = sizes
    else:
        sizes = blockSizes
    
    if type(blockStarts) is str:
        starts = [int(i) for i in blockStarts.split(',')]
        blockStarts = starts
    else:
        starts = blockStarts
    
    s = int(chromStart)
    block_ranges = [(s+start, s+start+size) for start,size in zip(starts,sizes)]
    return block_ranges

def explode_block_ranges(block_ranges):
    """Converts a list of block ranges generated by get_block_ranges()
    back into a set of three variables: chromStart, blockStarts, blockSizes"""
    chromStart = block_ranges[0][0] # Leftmost position
    blockStarts = [a-chromStart for a,b in block_ranges]
    blockSizes = [b-a+1 for a,b in block_ranges]
    return chromStart, blockStarts, blockSizes

def parse_ELR_line(elr_line):
    """Parses one line of an ELR file into an ELdata namedtuple.
    Example:
      0 6787 - E282A87D294A113D86A586D90A91D48A129D144S 0 1.0
    """
    chrom_num, chromStart, strand, EL_CIGAR, source_num, weight_string = elr_line.rstrip().split('\t')
    chrom = int(chrom_num)
    source = int(source_num)
    weight = float(weight_string)
    chromStart = int(chromStart)
    label_indices = [i for i,character in enumerate(EL_CIGAR) if not character.isdigit()]
    feature_lengths = [int(EL_CIGAR[a+1:b]) for a,b in zip(label_indices[:-1],label_indices[1:])]
    ranges = []
    splice = []
    gap = False
    position = chromStart
    for i,f in enumerate(feature_lengths):
        if not gap:
            rightside = position + f
            ranges += [(position, rightside)]
            position = rightside
        else:
            position += f
            if EL_CIGAR[label_indices[i]] in ['A','D']:
                splice += [True]
            else:
                splice += [False]
        
        gap = not gap
    
    s_tag = e_tag = capped = False
    if EL_CIGAR[0] == 'C':
        s_tag = capped = True
    elif EL_CIGAR[0] == 'S':
        s_tag = True
    elif EL_CIGAR[-1] == 'C':
        s_tag = capped = True
    elif EL_CIGAR[-1] == 'S':
        s_tag = True
    
    if EL_CIGAR[0] == 'E':
        e_tag = True
    elif EL_CIGAR[-1] == 'E':
        e_tag = True
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)

def parse_BED_line(bed_line, chrom_dict, source_dict, source_string=None, s_tag=False, e_tag=False, capped=False, gaps_are_junctions=False):
    """Parses one line of a 12- or 15-column BED file into an ELdata namedtuple.
    Examples:
      Ath_chr1 6787 8737 AT1G01020.6 . - 6787 8737 0,0,0 6 282,294,86,90,48,144 0,369,776,1448,1629,1806
      Ath_chr1 6787 8737 . . - 0 0 109,74,116 6 282,294,86,90,48,144 0,369,776,1448,1629,1806 1.0 TAIR10.40 EADADADADADS
    """
    bed_elements = bed_line.rstrip().split('\t')
    label = None
    if len(bed_elements) == 15: # 15-column ELR-BED
        chrom_string, chromStart, end, readname, score, strand, mmnum, mmorder, rgb, blocknum, blockSizes, blockStarts, weight, source_string, label = bed_elements
    elif len(bed_elements) == 12: # Standard BED12
        chrom_string, chromStart, end, readname, score, strand, mmnum, mmorder, rgb, blocknum, blockSizes, blockStarts = bed_elements
        try:
            weight = float(score)
        except:
            weight = float(1)        
    else: # Standard BED
        chrom_string, chromStart, end, readname, score, strand = bed_elements[:6]
        try:
            weight = float(score)
        except:
            weight = float(1)
        
        blockStarts = [0]
        blockSizes = [int(end)-int(chromStart)]
    
    ranges = get_block_ranges(chromStart, blockStarts, blockSizes)
    if label is not None: # Determine what kind of end labels exist based on the label
        s_tag = e_tag = capped = False
        if strand == '+':
            if label[0] == 'C':
                s_tag = capped = True
            elif label[0] == 'S':
                s_tag = True
            
            if label[-1] == 'E':
                e_tag = True
            
            splice = [True if i=='D' else False for i in label[1:-1:2]]
        elif strand == '-':
            if label[-1] == 'C':
                s_tag = capped = True
            elif label[-1] == 'S':
                s_tag = True
            
            if label[0] == 'E':
                e_tag = True
            
            splice = [True if i=='A' else False for i in label[1:-1:2]]
        else:
            splice = [False]*(len(ranges)-1)
    else:
        if gaps_are_junctions:
            splice = [True]*(len(ranges)-1)
        else:
            splice = [False]*(len(ranges)-1)
    
    if chrom_string in chrom_dict:
        chrom = chrom_dict[chrom_string]
    else:
        chrom = chrom_string
    
    if source_string is not None:
        if source_string in source_dict:
            source = source_dict[source_string]
        else:
            source = source_string
    else:
        source = None
    
    # chrom source strand ranges splice s_tag e_tag capped weight
    return ELdata(chrom, source, strand, ranges, splice, s_tag, e_tag, capped, weight)

def elr_to_readobject(elr_line):
    """Converts an ELR line to an RNAseqMapping object
    """
    input_data = parse_ELR_line(elr_line)
    output_object = RNAseqMapping(input_data)
    return output_object

##############################################
def which(x,value=True):
    return [a for a,b in enumerate(x) if b==value]

def order(list,reverse=False):
    """Returns a list of positions in list that contain
    the sorted values, without modifying the list.
    Behavior matches that of the R builtin 'order()'
    """
    list_order = [k for v,k in sorted([(v,k) for k,v in enumerate(list)],reverse=reverse)]
    return list_order

def table(x):
    """Returns a list of (item, abundance) doubles for a list of items"""
    y = list(set(x))
    z = [(i, sum([1 for j in x if j==i])) for i in y]
    return dict(z)

def combinations(items):
    """Returns a list of all k-combinations of items
    (credit @ninjagecko, StackOverflow)
    """
    if len(items) == 0:
        return []
    return (set(IT.compress(items,mask)) for mask in IT.product(*[[0,1]]*len(items)))

def sort_edges(list_of_edges):
    return [
        list_of_edges[k]
        for v,k in sorted([
            (v,k)
            for k,v in enumerate([
                int(a)
                for a,b,c in [
                    i.split(',')
                    for i in list_of_edges
                ]
            ])
        ])
    ]

def bridge_head(bridge,strand):
    """Returns the 5' end of the bridge"""
    if strand in ['+','>']:
        return int(bridge.split(',')[0])
    else:
        return int(bridge.split(',')[-2])

def bridge_tail(bridge,strand):
    """Returns the 3' end of the bridge"""
    if strand in ['+','>']:
        return int(bridge.split(',')[-2])
    else:
        return int(bridge.split(',')[0])

def bridges_are_compatible(bridge_a, bridge_b):
    """Returns a boolean, True if bridge_a and bridge_b
    are both overlapping and non-branching.
    """
    BA = bridge_a.split('|')
    BB = bridge_b.split('|')
    if len(set(BA+BB)) < len(BA)+len(BB):
        # The two bridges are partially overlapping
        branched = False
        JA = [i for i in BA if 'J' in i]
        for j in JA:
            jl,jr,jt = j.split(',')
            for a,b,c in [i.split(',') for i in BB]:
                if a == jl and c == 'U':
                    return False
                elif b==jr and c=='U':
                    return False
                elif a==jl and b!=jr and c=='J':
                    return False
                elif a!=jl and b==jr and c=='J':
                    return False

        JB = [i for i in BB if 'J' in i]
        for j in JB:
            jl,jr,jt = j.split(',')
            for a,b,c in [i.split(',') for i in BA]:
                if a == jl and c == 'U':
                    return False
                elif b==jr and c=='U':
                    return False
                elif a==jl and b!=jr and c=='J':
                    return False
                elif a!=jl and b==jr and c=='J':
                    return False

        return True
    else:
        return False

def overlap_groups(list_of_bridges):
    """ Identifies all groups of overlapping transcripts.
    Returns an array of tuples."""
    # start an empty array for the locus
    overlap_groups = []
    if len(list_of_bridges) == 0:
        return overlap_groups

    # sort the list of bridges primarily by increasing left position,
    # secondarily by decreasing right position.
    start_end_list = [
        (int(i.split(',')[0]),int(i.split(',')[-2]),i)
        for i in list_of_bridges
    ]
    start_end_list = sorted(
        sorted(
            start_end_list,key=lambda x:x[1],reverse=True
        ),key=lambda x:x[0],reverse=False
    )
    lpos = 0
    rpos = 0
    current_group = []
    for new_lpos,new_rpos,bridge in sorted(start_end_list):
        if rpos == 0:
            # initialize with the first bridge
            lpos = new_lpos
            rpos = new_rpos

        if new_lpos > rpos:
            # Dump the overlap group and start again
            overlap_groups.append(current_group)
            current_group = []
            lpos = new_lpos
            rpos = new_rpos
        elif new_rpos > rpos:
            rpos = new_rpos

        current_group.append(bridge)

    overlap_groups.append(current_group)
    return overlap_groups


def subtract_range(list_of_ranges,new_range):
    """Subtracts a range (L,R) from an existing
    list of ranges.
    """
    if new_range[1] - new_range[0] <= 0:
        return list_of_ranges

    updated_ranges = []
    for r in list_of_ranges:
        overlap = overlap_type(new_range,r)
        so = sum(overlap)
        if so > 0 and so < 4:
            # new_range overlaps with r
            if so == 1:
                # new_range hangs over the left edge of r
                if new_range[1] < r[1]:
                    updated_ranges += [(new_range[1],r[1])]
            elif so == 3:
                # new_range hangs over the right edge of r
                if new_range[0] > r[0]:
                    updated_ranges += [(r[0],new_range[0])]
            else:
                if overlap == [0,0,1,1]:
                    # new_range contains r
                    pass
                else:
                    # new_range is contained within r
                    if new_range[0] > r[0]:
                        updated_ranges += [(r[0],new_range[0]),(new_range[1],r[1])]
                    else:
                        updated_ranges += [(new_range[1],r[1])]
        else:
            updated_ranges += [r]

    return updated_ranges

def spans(pos,left,right):
    if left < pos and right > pos:
        return True
    else:
        return False


class DaisyChain():
    def __init__(self,chrom='', start=0, node_locs=[], node_types=[], node_scores=[]):
        """Initialize a DaisyChain object with an immutable list of
        nodes that make up the network. Each node should be given
        a type (options = S>,S<,E>,E<,D>,D<,A>,A<) and
        a score (float of the number of reads supporting this node)
        """
        self.chrom = chrom
        self.start = start
        self.total_reads = 0
        self.reads_assigned = 0
        self.node_locs = node_locs
        self.node_from_pos = {}
        for i in range(len(self.node_locs)):
            self.node_from_pos[self.node_locs[i]] = i
        
        self.number_of_nodes = len(self.node_locs)
        self.node_types = node_types
        self.node_scores = node_scores
        nodes_by_type = {}
        for node_type in ['S>','S<','E>','E<','D>','D<','A>','A<']:
            nodes_by_type[node_type] = [
                i for i in range(self.number_of_nodes)
                if self.node_types[i] == node_type
            ]
        
        self.nodes_by_type = nodes_by_type
        self.end_scores = {}
        for node_type in ['S>','S<','E>','E<']:
            self.end_scores[node_type] = sum([self.node_scores[i] for i in self.nodes_by_type[node_type]])
        
        self.edges = {}
        self.bridges = {}
        self.bridgeL = {}
        self.bridgeR = {}
        self.paths = {'+':{},'-':{}}
        self.transcripts = {}
    
    def __repr__(self):
        if self.number_of_nodes == 0:
            return '<DaisyChain object (empty)>'
        elif self.is_a_transcript():
            return '<DaisyChain transcript>'
        else:
            return '<DaisyChain object>'
    
    def display_edges(self):
        out_lines = []
        for k,v in self.edges.items():
            a,b,c = k.split(',')
            out_lines += ['{}.{} {}.{} {} {}'.format(
                a,self.node_types[int(a)],
                b,self.node_types[int(b)],
                c,[round(i,1) for i in v]
            )]
        
        return '\n'.join(out_lines)

    
    def render(self,bridge='',fit=None,color=True,window_width=64):
        """Renders a graphical representation of the DaisyChain object.
        """
        ground_path = '--'
        junction_start = ' \\'
        junction_end = '/ '
        junction_pass = '__'
        node_sym = {
            'S>':'\x1b[7;34;40mS>\x1b[0m',
            'S<':'\x1b[7;34;40mS<\x1b[0m',
            'E>':'\x1b[6;30;41mE>\x1b[0m',
            'E<':'\x1b[6;30;41mE<\x1b[0m',
        }
        def highlight(string,highlight_color='2;30;47'):
            return '\x1b[{}m{}\x1b[0m'.format(highlight_color,string)
        
        def space_is_empty(layer,pos1,pos2):
            return all([i=='  ' for i in rendering[layer][pos1:(pos2+1)]])
        
        def display_number(num):
            n = str(num)[-2:]
            if len(n) == 1:
                n = ' '+n
            
            return n
        
        if fit:
            node_locs = fit.node_locs
            node_types = fit.node_types
            nodes = range(len(node_locs))
            subnodes = [i for i in range(len(self.node_locs)) if self.node_locs[i] in node_locs]
        else:
            node_locs = self.node_locs
            node_types = self.node_types
            nodes = range(len(node_locs))
            subnodes = nodes
        
        width = len(nodes)*2 - 1
        layers = max(
            [len(self.neighbors(i,True,junction_flow='both')) for i in nodes] +
            [len(self.neighbors(i,False,junction_flow='both')) for i in nodes]
        ) + 1
        rendering = [['  ']*width for i in range(layers)]
        
        for i in subnodes:           
            if color and bridge:
                if str(i) == bridge.split(',')[0]:
                    rendering[0][i*2] = highlight(node_types[i])
                elif '|{},'.format(i) in bridge or ',{},'.format(i) in bridge:
                    rendering[0][i*2] = highlight(node_types[i])
                else:
                    rendering[0][i*2] = node_types[i]
            else:
                rendering[0][i*2] = node_types[i]
                
            if i+1 < len(nodes):
                ground_edge = '{},{},U'.format(i,i+1,'U')
                if ground_edge in self.edges:
                    if bridge and ground_edge in bridge:
                        rendering[0][(i*2)+1] = highlight(ground_path)
                    else:
                        rendering[0][(i*2)+1] = ground_path
            
            spliced_edges = self.neighbors(i,type='J',junction_flow='both')
            for j in spliced_edges:
                donor_acceptor_pair = (i*2,j*2)
                leftpos = min(donor_acceptor_pair)
                rightpos = max(donor_acceptor_pair)
                layer = 1
                junction_added = False
                junction_edge = '{},{},J'.format(min(i,j),max(i,j))
                while not junction_added:
                    if space_is_empty(layer,leftpos,rightpos):
                        if bridge and junction_edge in bridge:
                            rendering[layer][leftpos] = highlight(junction_start)
                            rendering[layer][rightpos] = highlight(junction_end)
                            for k in range(leftpos+1,rightpos):
                                rendering[layer][k] = highlight(junction_pass)
                        else:
                            rendering[layer][leftpos] = junction_start
                            rendering[layer][rightpos] = junction_end
                            for k in range(leftpos+1,rightpos):
                                rendering[layer][k] = junction_pass
                        
                        junction_added = True
                    else:
                        layer += 1
                        if layer >= len(rendering):
                            rendering += [['  ']*width]
        firstline = flatten([[display_number(i),'  '] for i in nodes])[:-1]
        all_lines = [firstline] + [
            j for j in rendering
            if not all([x == ' ' for x in j])
        ]
        
        lines_to_render = []
        ldivide = 0
        rdivide = window_width
        while ldivide < len(firstline):
            new_lines = [i[ldivide:rdivide] for i in all_lines]
            if color:
                new_lines[0] = '\x1b[6;37;40m'+''.join(new_lines[0])+'\x1b[0m'
                for k,v in node_sym.items():
                    for i in range(1,len(new_lines)):
                        new_lines[i] = ''.join(new_lines[i])
                        new_lines[i] = new_lines[i].replace(k,v)
            
            lines_to_render += ['\n'.join(new_lines)]
            ldivide += window_width
            rdivide += window_width
        
        if color:
            lines_to_render = ['\x1b[2;30;47m{}:{}-{}\x1b[0m'.format(self.chrom,self.start,self.start+self.node_locs[-1])]+lines_to_render
        else:
            lines_to_render = ['{}:{}-{}'.format(self.chrom,self.start,self.start+self.node_locs[-1])]+lines_to_render
        
        return '\n'.join(lines_to_render)
        
    
    def add_edge(self,left_node,right_node,type,weights):
        """Adds an item to the dict 'edges' that connects two nodes.
        Properties of an 'edges' item:
          key = 'left_node,right_node,type'
          value = [forward_weight,reverse_weight,nonstranded_weight]
        """
        self.paths = {'+':{},'-':{}}
        edge_key = '{},{},{}'.format(left_node, right_node, type)
        if edge_key in self.edges:
            current_weights = self.edges[edge_key]
            updated_weights = [a+b for a,b in zip(current_weights,weights)]
            updated_weights = [i if i >= 0 else 0 for i in updated_weights]
            if sum(updated_weights) > 0:
                self.edges[edge_key] = updated_weights
            else:
                self.remove_edge(edge_key)
        else:
            if all([i>=0 for i in weights]) and sum(weights) > 0:
                self.edges[edge_key] = weights
    
    def remove_edge(self,edge_key):
        """Deletes an edge from the 'edges' dict if it exists.
        """
        self.paths = {'+':{},'-':{}}
        if edge_key in self.edges:
            del self.edges[edge_key]
            for bridge in self.bridges.keys():
                if edge_key in bridge:
                    self.bridges[bridge] = [0,0,0]
    
    def add_bridge(self,list_of_edge_keys,weights):
        bridge_key = '|'.join(list_of_edge_keys)
        if bridge_key in self.bridges:
            current_weights = self.bridges[bridge_key]
            updated_weights = [a+b for a,b in zip(current_weights,weights)]
            updated_weights = [i if i >= 0 else 0 for i in updated_weights]
            self.bridges[bridge_key] = updated_weights
        else:
            weights = [i if i >= 0 else 0 for i in weights]
            if sum(weights) > 0:
                self.bridges[bridge_key] = weights
                leftmost = int(list_of_edge_keys[0].split(',')[0])
                rightmost = int(list_of_edge_keys[-1].split(',')[1])
                if leftmost in self.bridgeL:
                    self.bridgeL[leftmost].add(bridge_key)
                else:
                    self.bridgeL[leftmost] = set([bridge_key])
                
                if rightmost in self.bridgeR:
                    self.bridgeR[rightmost].add(bridge_key)
                else:
                    self.bridgeR[rightmost] = set([bridge_key])
    
    def explode_edges(self, include_a=None, include_b=None):
        """Returns a list of triples for all edges in the network:
        (left_node, right_node, edge_type)
        Can require a specific node, two specific nodes, or no
        specific nodes to be in the edge.
        """
        edge_strings = [key.split(',') for key in sorted(list(self.edges.keys()))]
        exploded_edges = [(int(a),int(b),c) for a,b,c in edge_strings]
        if include_a is None and include_b is None:
            return exploded_edges
        else:
            if include_b is None:
                return [i for i in exploded_edges if include_a in i]
            elif include_a is None:
                return [i for i in exploded_edges if include_b in i]
            else:
                return [i for i in exploded_edges if include_a in i and include_b in i]

    def neighbors(self,node_index,forward_strand=True,junction_flow='downstream',type='all',return_type='node'):
        """Returns a list of node indices that can be
        reached via an edge from the provided
        node_index. If forward_strand, only searches
        to the right, else only to the left.
        By default, a node connected only by a J-type
        edge in the wrong orientation will not be returned.
        """
        connected_edges = self.explode_edges(node_index)
        node_type = self.node_types[node_index]
        if junction_flow is None or type=='U':
            connected_edges = [i for i in connected_edges if i[2] != 'J']
        
        if type == 'J':
            connected_edges = [i for i in connected_edges if i[2] != 'U']
        
        if forward_strand:
            if (junction_flow == 'downstream' and node_type != 'D>') \
            or (junction_flow == 'upstream' and node_type != 'A<'):
                connected_edges = [i for i in connected_edges if i[2] != 'J']
            
            if return_type == 'node':
                connected_nodes = [b for a,b,c in connected_edges if a == node_index]
            else:
                connected_edges = [i for i in connected_edges if i[0] == node_index]
        else:
            if (junction_flow == 'downstream' and node_type != 'D<') \
            or (junction_flow == 'upstream' and node_type != 'A>'):
                connected_edges = [i for i in connected_edges if i[2] != 'J']
            
            if return_type == 'node':
                connected_nodes = [a for a,b,c in connected_edges if b == node_index]
            else:
                connected_edges = [i for i in connected_edges if i[1] == node_index]
        
        if return_type == 'node':
            return connected_nodes
        else:
            return [','.join([str(j) for j in i]) for i in connected_edges]
    
    def valid_node_pairs_exist(self):
        """Boolean indicating whether it would be possible to assemble
        at least one transcript end-to-end from the available
        nodes. This would require an ordered S>E> pair or an E<S< pair.
        """
        S_plus = False
        E_plus = False
        S_minus = False
        E_minus = False
        for node_type in self.node_types:
            if node_type == 'S>':
                S_plus = True
            if node_type == 'E<':
                E_minus = True
            if node_type == 'E>' and S_plus:
                return True
            if node_type == 'S<' and E_minus:
                return True
        
        return False
    
    def path_exists(self,source,sink,forward_strand=None):
        """Returns Boolean indicating whether at least one valid 
        path exists from source to sink. Updates the dict
        self.paths, so that if this source->sink combination
        had already been checked before, it will not check twice.
        """
        junction_flow='downstream'
        if forward_strand is None:
            if sink < source:
                forward_strand = False
                strand = '-'
            else:
                forward_strand = True
                strand = '+'
        else:
            if forward_strand == True:
                strand = '+'
            else:
                strand = '-'
        
        if forward_strand and source > sink:
            return False
        elif not forward_strand and source < sink:
            return False
        
        if source not in self.paths[strand]:
            self.paths[strand][source] = {}
        
        if sink in self.paths[strand][source]:
            return self.paths[strand][source][sink]
        
        if source == sink:
            self.paths[strand][source][sink] = True
            return True
        
        connections = self.neighbors(source,forward_strand,junction_flow)
        if sink in connections:
            self.paths[strand][source][sink] = True
            return True
        else:
            for downstream_node in connections:
                if (forward_strand and downstream_node < sink) or \
                (not forward_strand and downstream_node > sink):
                    if self.path_exists(downstream_node,sink,forward_strand):
                        self.paths[strand][source][sink] = True
                        return True
        
        self.paths[strand][source][sink] = False
        return False
    
    def path_is_unique(self,source,sink):
        """Returns Boolean indicating whether exactly
        one path connects source to sink.
        """
        if source == sink:
            return True
        
        if sink < source:
            forward_strand = False
            strand = '-'
        else:
            forward_strand = True
            strand = '+'
        
        active_node = source
        while active_node != sink:
            connections = self.neighbors(active_node,forward_strand)
            if len(connections) == 1:
                active_node = connections[0]
            else:
                paths_exist = [self.path_exists(i,sink,forward_strand) for i in connections]
                if sum(paths_exist) > 1:
                    return False
                else:
                    active_node = connections[which(paths_exist)[0]]
        
        return True
    
    def contains_paths(self):
        """ Returns True if at least one path
        exists between at least one S node
        and at least one same-stranded E node.
        """
        S_plus = self.nodes_by_type['S>']
        S_minus = self.nodes_by_type['S<']
        E_plus = self.nodes_by_type['E>']
        E_minus = self.nodes_by_type['E<']
        
        # Generate all pairs of nodes on the same strand in a compatible direction
        putative_pairs_plus = flatten([[(i,j) for j in E_plus if j > i] for i in S_plus])
        putative_pairs_minus = flatten([[(i,j) for j in E_minus if j < i] for i in S_minus])
        if len(putative_pairs_plus) == 0 and len(putative_pairs_minus) == 0:
            return False
        
        for p in putative_pairs_plus + putative_pairs_minus:
            if self.path_exists(p[0],p[1]):
                return True
        
        return False
    
    def is_a_transcript(self):
        """Boolean indicating whether the DaisyChain object satisfies:
          1) The nodes on either end are a compatible pair of start/end nodes
          2) One and only one path exists from the start to end nodes
        If true, the object is given four new attributes:
          transcript_length = number of nucleotides in the transcript model
          transcript_strand = +,-
          transcript_S = genomic coordinates of transcript start site
          transcript_E = genomic coordinates of transcript end site
        """
        first_node_type = self.node_types[0]
        last_node_type = self.node_types[-1]
        if (first_node_type,last_node_type) == ('S>','E>'):
            start_node = 0
            end_node = self.number_of_nodes - 1
            forward_strand = True
        elif (first_node_type,last_node_type) == ('E<','S<'):
            start_node = self.number_of_nodes - 1
            end_node = 0
            forward_strand = False
        else:
            return False
        
        current_node = start_node
        while current_node != end_node:
            downstream_nodes = self.neighbors(current_node,forward_strand)
            if len(downstream_nodes) == 1:
                current_node = downstream_nodes[0]
            else:
                return False
        
        exploded_edges = self.explode_edges()
        self.transcript_length = sum([
            self.node_locs[b]-self.node_locs[a]
            for a,b,c in exploded_edges if c == 'U'
        ])
        if forward_strand:
            self.transcript_strand =  '+'
            self.transcript_S = self.node_locs[0] + self.start
            self.transcript_E = self.node_locs[-1] + self.start
        else:
            self.transcript_strand =  '-'
            self.transcript_S = self.node_locs[-1] + self.start
            self.transcript_E = self.node_locs[0] + self.start
        
        return True
    
    def stranded_ratio(self,node):
        """Returns a float of the proportion of flow
        estimated to be plus stranded. The total score
        of all downstream sinks and upstream sources are
        summed for both strands, producing a ratio between
        0 (100% minus) and 1 (100% plus)
        """
        plus_flow = 0
        minus_flow = 0

        total_S = self.end_scores['S>'] + self.end_scores['S<']
        total_E = self.end_scores['E>'] + self.end_scores['E<']

        up_S_plus = sum([self.node_scores[i] for i in self.nodes_by_type['S>'] if self.path_exists(i,node,True)])/total_S
        up_S_minus = sum([self.node_scores[i] for i in self.nodes_by_type['S<'] if self.path_exists(i,node,False)])/total_S
        down_E_plus = sum([self.node_scores[i] for i in self.nodes_by_type['E>'] if self.path_exists(node,i,True)])/total_E
        down_E_minus = sum([self.node_scores[i] for i in self.nodes_by_type['E<'] if self.path_exists(node,i,False)])/total_E

        total_stranded = sum([up_S_plus,up_S_minus,down_E_plus,down_E_minus])
        total_plus = up_S_plus + down_E_plus
        if total_stranded == 0:
            return 0.5
        
        ratio = total_plus/total_stranded
        return ratio
    
    def sinks_from_source(self,source):
        """Returns the position of all sink nodes
        (E>,E<) reachable from a source node.
        Includes the source node if source is a sink.
        """
        source_type = self.node_types[source]
        sinks = []
        if source_type in ['E>','E<']:
            sinks = [source]
        
        if '>' in source_type:
            sinks += [i for i in self.nodes_by_type['E>'] if i > source]
        elif '<' in source_type:
            sinks += [i for i in self.nodes_by_type['E<'] if i < source]
        
        return [i for i in sinks if self.node_scores[i] > 0 and self.path_exists(source,i)]
    
    def sources_from_sink(self,sink):
        """Returns the position of all source nodes
        (S>,S<) reachable from a sink node, 
        with junction flow reversed.
        """
        sink_type = self.node_types[sink]
        sources = []
        if sink_type in ['S>','S<']:
            sources = [sink]
        
        if '>' in sink_type:
            sources += [i for i in self.nodes_by_type['S>'] if i < sink]
        elif '<' in sink_type:
            sources += [i for i in self.nodes_by_type['S<'] if i > sink]
        
        return [i for i in sources if self.node_scores[i] > 0 and self.path_exists(i,sink)]
    
    def proportions(self,node,length_normalize=False,return_weight=False,compatible=''):
        """Takes any node and defines the proportion of
        reads contained by edges from that node.
        For S/E nodes, all bridges upstream vs. downstream are compared.
        For D/A nodes, alternate paths are different J edges or U edges
        flowing downstream from D or upstream from A.
        Returns a dict of {'edge_ID':proportion,...} that sums to 1.
        """
        edge_weights = {}
        total = 0
        node_type = self.node_types[node]
        relevant_bridges = set(
            [i for i in self.bridges.keys() if ',{},'.format(node) in i] \
            + list(self.bridgeL.get(node,[])) \
            + list(self.bridgeR.get(node,[]))
        )
        
        relevant_bridges = [i for i in relevant_bridges if i in self.bridges]
        if compatible:
            # Trim relevant_bridges to only those that
            # are compatible with the given bridge
            relevant_bridges = [i for i in relevant_bridges if bridges_are_compatible(i,compatible)]
        
        if node_type[0] in ['S','E']:
            # For S/E nodes, .proportions() reports the estimated
            # proportion of reads 'captured' by the node vs. flowthrough
            if node_type in ['S>','E>']:
                upstream_edge = '{},{},U'.format(node-1,node)
                downstream_edge = '{},{},U'.format(node,node+1)
            else:
                downstream_edge = '{},{},U'.format(node-1,node)
                upstream_edge = '{},{},U'.format(node,node+1)
            
            if node_type in ['S>','E<']:
                upstream_bridges = [i for i in relevant_bridges if upstream_edge in i]
                downstream_bridges = [i for i in relevant_bridges if downstream_edge in i and i not in upstream_bridges]
            else:
                downstream_bridges = [i for i in relevant_bridges if downstream_edge in i]
                upstream_bridges = [i for i in relevant_bridges if upstream_edge in i and i not in downstream_bridges]
            
            if length_normalize:
                upstream_bridge_lengths = [float(self.bridge_length(i)) for i in upstream_bridges]
                downstream_bridge_lengths = [float(self.bridge_length(i)) for i in downstream_bridges]
                upstream_weights = [
                    [v[0]/l,v[1]/l,v[2]/l]
                    for v,l in zip(
                        [self.bridges[i] for i in upstream_bridges],
                        upstream_bridge_lengths
                    ) if l > 0
                ]
                downstream_weights = [
                    [v[0]/l,v[1]/l,v[2]/l]
                    for v,l in zip(
                        [self.bridges[i] for i in downstream_bridges],
                        downstream_bridge_lengths
                    ) if l > 0
                ]
            else:
                upstream_weights = [self.bridges[i] for i in upstream_bridges]
                downstream_weights = [self.bridges[i] for i in downstream_bridges]
            
            strand_ratio = self.stranded_ratio(node)
            
            total = sum(flatten(upstream_weights))+sum(flatten(downstream_weights))
            if node_type == 'S>':
                captured = sum([p+(strand_ratio*n) for p,m,n in downstream_weights])
            elif node_type == 'E>':
                captured = sum([p+(strand_ratio*n) for p,m,n in upstream_weights])
            elif node_type == 'S<':
                captured = sum([m+((1-strand_ratio)*n) for p,m,n in downstream_weights])
            elif node_type == 'E<':
                captured = sum([m+((1-strand_ratio)*n) for p,m,n in upstream_weights])
            
            passed = total - captured            
            edge_weights = {'captured':captured,'passed':passed}
        
        elif node_type == 'D>':
            out_edges = self.neighbors(node,forward_strand=True,junction_flow='downstream',return_type='edge')
            for edge in out_edges:
                in_edge = [i for i in relevant_bridges if edge in i]
                if length_normalize:
                    in_edge_lengths = [float(self.bridge_length(i)) for i in in_edge]
                    in_weights = [
                        [v[0]/l,v[1]/l,v[2]/l]
                        for v,l in zip(
                            [self.bridges[i] for i in in_edge],
                            in_edge_lengths
                        ) if l > 0 
                    ]
                else:
                    in_weights = [self.bridges[i] for i in in_edge]
                
                total += sum(flatten(in_weights))
                edge_weights[edge] = sum([p+n for p,m,n in in_weights])
        elif node_type == 'D<':
            out_edges = self.neighbors(node,forward_strand=False,junction_flow='downstream',return_type='edge')
            for edge in out_edges:
                in_edge = [i for i in relevant_bridges if edge in i]
                if length_normalize:
                    in_edge_lengths = [float(self.bridge_length(i)) for i in in_edge]
                    in_weights = [
                        [v[0]/l,v[1]/l,v[2]/l]
                        for v,l in zip(
                            [self.bridges[i] for i in in_edge],
                            in_edge_lengths
                        ) if l > 0
                    ]
                else:
                    in_weights = [self.bridges[i] for i in in_edge]
                
                total += sum(flatten(in_weights))
                edge_weights[edge] = sum([m+n for p,m,n in in_weights])
        elif node_type == 'A>':
            in_edges = self.neighbors(node,forward_strand=False,junction_flow='upstream',return_type='edge')
            for edge in in_edges:
                out_edge = [i for i in relevant_bridges if edge in i]
                if length_normalize:
                    out_edge_lengths = [float(self.bridge_length(i)) for i in out_edge]
                    out_weights = [
                        [v[0]/l,v[1]/l,v[2]/l]
                        for v,l in zip(
                            [self.bridges[i] for i in out_edge],
                            out_edge_lengths
                        ) if l > 0
                    ]
                else:
                    out_weights = [self.bridges[i] for i in out_edge]
                
                total += sum(flatten(out_weights))
                edge_weights[edge] = sum([p+n for p,m,n in out_weights])
        elif node_type == 'A<':
            in_edges = self.neighbors(node,forward_strand=True,junction_flow='upstream',return_type='edge')
            for edge in in_edges:
                out_edge = [i for i in relevant_bridges if edge in i]
                if length_normalize:
                    out_edge_lengths = [float(self.bridge_length(i)) for i in out_edge]
                    out_weights = [
                        [v[0]/l,v[1]/l,v[2]/l]
                        for v,l in zip(
                            [self.bridges[i] for i in out_edge],
                            out_edge_lengths
                        ) if l > 0
                    ]
                else:
                    out_weights = [self.bridges[i] for i in out_edge]
                
                total += sum(flatten(out_weights))
                edge_weights[edge] = sum([m+n for p,m,n in out_weights])
        
        edge_weights = dict([(k,v) for k,v in edge_weights.items() if v > 0])
        
        if return_weight:
            return edge_weights
        
        if total == 0:
            return dict()
        
        return dict([(k,v/total) for k,v in edge_weights.items()])
        
    
    def prune_edges(self,proportion=0.05):
        """Removes any edges that capture less than 'proportion'
        of the reads flowing ending at a given node.
        """
        for i in range(0,self.number_of_nodes):
            node_type = self.node_types[i]
            edge_proportions = self.proportions(i)
            # print(i,node_type,edge_proportions)
            if node_type[0] in ['S','E']:
                if edge_proportions.get('passed',0) < proportion:
                    if node_type in ['S>','E>']:
                        upstream_edge = '{},{},U'.format(i-1,i)
                        downstream_edge = '{},{},U'.format(i,i+1)
                    else:
                        upstream_edge = '{},{},U'.format(i,i+1)
                        downstream_edge = '{},{},U'.format(i-1,i)
                    
                    if node_type in ['S>','S<']:
                        self.remove_edge(upstream_edge)
                    else:
                        self.remove_edge(downstream_edge)
            else:
                for k,v in edge_proportions.items():
                    if v < proportion:
                        self.remove_edge(k)
    
    def flowthrough(self,node):
        """For a source or a sink feature, flowthrough is the total score
        of the given source/sink node + all source/sink nodes beyond the
        given node to which a path exists.
        """
        node_type = self.node_types[node]
        if node_type not in ['S>','S<','E>','E<']:
            return 0
        node_score = self.node_scores[node]
        if 'S' in  node_type:
            upstream_sources = self.sources_from_sink(node)[1:]
            source_scores = sum([node_score]+[self.node_scores[i] for i in upstream_sources])
            return source_scores
        elif 'E' in node_type:
            downstream_sinks = self.sinks_from_source(node)[1:]
            sink_scores = sum([node_score]+[self.node_scores[i] for i in downstream_sinks])
            return sink_scores
        
    def dominant_pair(self):
        """Returns a path-connected pair of ends that...
        (1) has the lowest flowthrough, and
        (2) has the dominant signal
        """
        best_score = 0
        dominant_pair = None
        for i in self.nodes_by_type['S>'] + self.nodes_by_type['S<']:
            source_score = float(self.node_scores[i])
            if source_score == 0:
                continue
            
            flowthrough = self.proportions(i).get('captured',0)
            flowsource = source_score * flowthrough
            sinks = sorted(self.sinks_from_source(i))
            if len(sinks) == 0:
                continue
            
            for j in sinks:
                sink_score = float(self.node_scores[j])
                if sink_score == 0:
                    continue
                
                flowthrough = self.proportions(j).get('captured',0)
                flowsink = sink_score * flowthrough
                test_score = flowsource + flowsink
                if test_score > best_score:
                    # Update the best score and dominant pair
                    dominant_pair = (i,j)
                    best_score = test_score
                elif test_score == best_score and test_score > 0:
                    # Evaluate if the current pair is farther apart
                    # than the previous pair
                    if abs(self.node_locs[i] - self.node_locs[j]) > \
                    abs(self.node_locs[dominant_pair[0]] - self.node_locs[dominant_pair[1]]):
                        dominant_pair = (i,j)
                        best_score = test_score
        
        return dominant_pair

    def heaviest_path(self,source,sink,strand):
        """Finds the heaviest supported
        non-branching path through edges
        from source to sink.
        """
        if not self.path_exists(source,sink):
            return
        
        strand_ratios = [round(self.stranded_ratio(i),2) for i in range(self.number_of_nodes)]
        if strand == '+':
            if source > sink:
                return
            forward_strand = True
        elif strand == '-':
            if source < sink:
                return
            forward_strand = False
        
        source_prop = self.proportions(source).get('captured',1.0)
        sink_prop = self.proportions(sink).get('captured',1.0)
        min_prop = min([source_prop,sink_prop])
        
        heaviest_path = []
        active_node = source
        while active_node != sink:
            neighbors = self.neighbors(active_node,forward_strand,junction_flow='downstream',type='all',return_type='edge')
            
            if len(neighbors) == 1:
                heaviest_path.append(neighbors[0])
                l,r,t = neighbors[0].split(',')
                if strand == '+':
                    active_node = int(r)
                elif strand == '-':
                    active_node = int(l)
            else:
                props = sorted([(v,k) for k,v in self.proportions(active_node,compatible='|'.join(heaviest_path)).items()],reverse=True)
                try:
                    props = [(v,k) for v,k in props if self.path_exists(bridge_tail(k,strand),sink,forward_strand)]
                except:
                    print('ERROR: heaviest_path() failed.')
                    print(props)
                    print(self.render())
                    print(source,sink,active_node,neighbors)
                    sys.exit(1)
                
                if len(props) == 0:
                    return
                
                for v,k in props:
                    if k in neighbors:
                        break
                
                heaviest_path.append(k)
                l,r,t = k.split(',')
                if strand == '+':
                    active_node = int(r)
                elif strand == '-':
                    active_node = int(l)                
                
                if v < min_prop:
                    min_prop = v
        
        heaviest_path = sort_edges(heaviest_path)
        return ('|'.join(heaviest_path),strand,min_prop)
    
    def bridges_in_range(self,left_node,right_node,overflow=True):
        """Returns a set of all bridges that are contained
        in the space between two nodes
        """
        l = int(left_node)
        r = int(right_node)
        bridges_in_range = [
            i for i in self.bridges.keys()
            if int(i.split(',')[0]) >= l
            and int(i.split(',')[-2]) <= r
        ]
        if overflow:
            bridges_in_range += [
                i for i in self.bridges.keys()
                if ',{},'.format(l) in i
                or ',{},'.format(r) in i
            ]
            bridges_in_range = list(set(bridges_in_range))
        
        return bridges_in_range
    
    def contained_bridges(self,bridge):
        """Takes a bridge string and returns
        a list of bridges in the network that
        are subsets of the input bridge.
        """
        leftmost = bridge.split(',')[0]
        rightmost = bridge.split(',')[-2]
        bridges_in_range = self.bridges_in_range(leftmost,rightmost,overflow=False)
        contained_bridges = [
            i for i in bridges_in_range
            if bridges_are_compatible(bridge,i)
        ]
        return contained_bridges
    
    def bridge_length(self,bridge):
        """Returns the summed distance in nucleotides
        between all nodes for U edges in a bridge.
        """
        ground_edges = [i.split(',') for i in bridge.split('|') if 'U' in i]
        lengths = [1 + abs(self.node_locs[int(b)] - self.node_locs[int(a)]) for a,b,c in ground_edges]
        return sum(lengths)
        
    
    def extract_path(self,path,strand,proportion):
        """Modifies the edges and bridges of the object
        to pull out 'proportion' of all reads from the
        paths in question. Returns the number of reads
        assigned to the extracted path.
        """
        # path_edge_weights = dict([(i,float(0)) for i in path.split('|')])
        reads_assigned = float(0)
        bridges_altered = []
        leftmost = int(path.split(',')[0])
        rightmost = int(path.split(',')[-2])
        if strand == '+':
            source = leftmost
            sink = rightmost
        else:
            source = rightmost
            sink = leftmost
        
        if self.path_is_unique(source,sink):
            other_sinks = False
            other_sources = False
            if len(self.sinks_from_source(source)) > 1:
                # Deplete proportion relative to sink flowthrough
                # proportion = proportion * (self.node_scores[sink] / self.flowthrough(sink))
                other_sinks = True
            
            if len(self.sources_from_sink(sink)) > 1:
                # proportion = proportion * (self.node_scores[source] / self.flowthrough(source))
                other_sources = True
            
            if not other_sinks and not other_sources:
                proportion = 1
                proportion_to_extract = 1
                self.node_scores[leftmost] = 0
                self.node_scores[rightmost] = 0
        
        total_spanning_reads = sum([sum(self.bridges[i]) for i in self.bridges_in_range(leftmost,rightmost,overflow=True)])
        if total_spanning_reads == 0:
            for edge in path.split('|'):
                self.remove_edge(edge)
            
            return (0,[])
        
        contained_bridges = self.contained_bridges(path)
        for bridge in contained_bridges:
            if proportion != 1:
                l = int(bridge.split(',')[0])
                r = int(bridge.split(',')[-2])
                spanning_reads = sum([sum(self.bridges[i]) for i in self.bridges_in_range(l,r,overflow=False)])
                if strand == '+':
                    stranded_compatible_reads = sum([self.bridges[i][0]+self.bridges[i][2] for i in self.contained_bridges(bridge)])
                else:
                    stranded_compatible_reads = sum([self.bridges[i][1]+self.bridges[i][2] for i in self.contained_bridges(bridge)])
                
                if stranded_compatible_reads == 0 or spanning_reads == 0:
                    proportion_to_extract = 1
                else:
                    compatible_proportion = stranded_compatible_reads / spanning_reads
                    proportion_to_extract = compatible_proportion / proportion
                
                if proportion_to_extract > 1:
                    proportion_to_extract = 1
            else:
                proportion_to_extract = 1
            
            negative_weight_triple = [-float(i)*proportion_to_extract for i in self.bridges[bridge]]
            if strand == '+':
                negative_weight_triple[1] = 0
            elif strand == '-':
                negative_weight_triple[0] = 0
            
            weight = -sum(negative_weight_triple)
            # Add weights to path and remove from network
            reads_assigned += weight
            split_bridge = bridge.split('|')
            self.add_bridge(split_bridge,negative_weight_triple)
            bridges_altered += [bridge]
            for edge in split_bridge:
                l,r,t = edge.split(',')
                self.add_edge(l,r,t,negative_weight_triple)
                # path_edge_weights[edge] += weight
        
        proportion_assigned = reads_assigned / total_spanning_reads
        if reads_assigned == 0:
            self.node_scores[source] = 0
            self.node_scores[sink] = 0
        
        self.node_scores[leftmost] = self.node_scores[leftmost] - (self.node_scores[leftmost]*proportion)
        self.node_scores[rightmost] = self.node_scores[rightmost] - (self.node_scores[rightmost]*proportion)
        return (reads_assigned,bridges_altered)
    
    def maximum_flow(self,minimum_proportion=0.02,readlength=50):
        """Creates a minimized version of the DaisyChain object
        and iteratively extracts the heaviest path until no more
        transcripts can be assembled. Stores these transcript models
        in the DaisyChain object's 'transcripts' attribute.
        """
        ratios = [round(self.stranded_ratio(i),2) for i in range(self.number_of_nodes)]
        edge_triples = self.edges
        stranded_capacity = dict()
        for k,v in edge_triples.items():
            # Assign ns weights to p and m strands by ratios
            split_edge = k.split(',')
            l = int(split_edge[0])
            r = int(split_edge[1])
            edge_ratio = (ratios[l] + ratios[r])/2
            if split_edge[-1] == 'U':
                # readlength_bonus = readlength*2
                capacity_plus = (v[0] + edge_ratio*v[2])/(self.node_locs[r]-self.node_locs[l])
                capacity_minus = (v[1] + (1-edge_ratio)*v[2])/(self.node_locs[r]-self.node_locs[l])
            elif split_edge[-1] == 'J':
                # readlength_bonus = readlength
                capacity_plus = (v[0] + edge_ratio*v[2])
                capacity_minus = (v[1] + (1-edge_ratio)*v[2])
            
            stranded_capacity[k] = (capacity_plus, capacity_minus)
        
        # return(stranded_capacity)
        
        dc = copy.deepcopy(self)
        dc.transcripts = {}
        dc.paths = {'+':{},'-':{}}
        extracted_transcripts = {}
        ends_to_use = dc.dominant_pair()
        while ends_to_use:
            source, sink = ends_to_use
            # print(dc.render())
            # print(ends_to_use)
            if source is None or sink is None:
                ends_to_use = None
                continue
            
            if source < sink:
                strand = '+'
            elif source > sink:
                strand = '-'
            else:
                ends_to_use == None
                continue
            
            path = dc.heaviest_path(ends_to_use[0],ends_to_use[1],strand)
            if path:
                reads_assigned,bridges_altered = dc.extract_path(path[0],path[1],path[2])
                # print(reads_assigned)
                if reads_assigned:
                    edges_checked = set()
                    for bridge in bridges_altered:
                        if sum(dc.bridges[bridge]) < sum(self.bridges[bridge]) * minimum_proportion:
                            dc.bridges[bridge] = [0,0,0]
                        
                        for edge in bridge.split('|'):
                            if edge not in edges_checked:
                                if sum(dc.edges.get(edge,[])) < sum(self.edges.get(edge,[])) * minimum_proportion:
                                    dc.remove_edge(edge)
                                
                                edges_checked.add(edge)
                    
                    if path[0] in extracted_transcripts:
                        extracted_transcripts[path[0]]['reads_assigned'] += reads_assigned
                    else:
                        extracted_transcripts[path[0]] = {
                            'strand':path[1],
                            'reads_assigned':reads_assigned
                        }
            else:
                dc.node_scores[ends_to_use[0]] = 0
                dc.node_scores[ends_to_use[1]] = 0
            
            ends_to_use = dc.dominant_pair()
        
        self.transcripts = extracted_transcripts
    
    def EM_read_assignment(self,leftmost,rightmost,transcripts,priors,terminal_delta=1,digits=3):
        """Given a set of transcript models and prior
        assigned weights, iteratively reassigns bridge
        weights to the transcripts, updating priors until
        'delta' (total change in coverage between iterations)
        falls below the given cap.
        Returns two values:
            -A locus score, defined as 
            -A list of updated read assignments
        """
        # Identify all bridges in the relevant range
        relevant_bridges = set(self.bridges_in_range(leftmost,rightmost,overflow=True))
        total_reads = sum([sum(self.bridges[i]) for i in relevant_bridges])
        priors_dict = dict(zip(transcripts,priors))
        transcript_edges = {}
        for t in transcripts:
            transcript_edges[t] = {}
            for e in t.split('|'):
                transcript_edges[t][e] = float(0)
        
        delta = float(10)
        while delta > terminal_delta:
            # Repeat until the total change in read coverage between
            # iterations drops to terminal_delta
            previous_priors = copy.deepcopy(priors)
            assigned_reads = dict(zip(transcripts,[float(0)]*len(transcripts)))
            unassigned_bridges = dict([(i,copy.deepcopy(self.bridges[i])) for i in relevant_bridges])
            coverages = [priors_dict[i]/self.bridge_length(i) for i in transcripts]
            t_dict = {}
            # Find all bridges compatible with each transcript model
            for t in transcripts:
                t_dict[t] = set(self.contained_bridges(t))
            
            bridges = list(relevant_bridges)
            for bridge in bridges:
                if sum(unassigned_bridges[bridge]) == 0:
                    continue
                
                compatible_transcripts = [t for t in t_dict.keys() if bridge in t_dict[t]]
                if len(compatible_transcripts) == 0:
                    continue
                elif len(compatible_transcripts) == 1:
                    # Assign the full strand-compatible weight
                    # of the bridge to its compatible transcript
                    t = compatible_transcripts[0]
                    strand = self.transcripts[t]['strand']
                    weight_triple = unassigned_bridges[bridge]
                    if strand == '+':
                        weight_triple[1] = 0
                    elif strand == '-':
                        weight_triple[0] = 0
                    
                    reads_to_assign = sum(weight_triple)
                    t = compatible_transcripts[0]
                    assigned_reads[t] += reads_to_assign
                    for e in bridge.split('|'):
                        if e in transcript_edges[t]:
                            transcript_edges[t][e] += reads_to_assign
                    unassigned_bridges[bridge] = [a-b for a,b in zip(unassigned_bridges[bridge],weight_triple)]
                else:
                    # More than one compatible transcript model exits
                    # 1) Compute the proportion to assign each model based on previous priors
                    # 2) Assign the full weight of the bridge to all compatible models
                    weight_triple = unassigned_bridges[bridge]
                    strands = [self.transcripts[t]['strand'] for t in compatible_transcripts]
                    strand_number = table(strands)
                    
                    # Assign strand-specific reads to the subset of models 
                    if weight_triple[0]:
                        if strand_number.get('+',0) == 0:
                            weight_triple[0] = 0
                        elif strand_number.get('+',0) == 1:
                            reads_to_assign = weight_triple[0]
                            t = compatible_transcripts[which(strands,'+')[0]]
                            assigned_reads[t] += reads_to_assign
                            for e in bridge.split('|'):
                                if e in transcript_edges[t]:
                                    transcript_edges[t][e] += reads_to_assign
                        else:
                            stranded_subset = [(t,c) for t,c,s in zip(compatible_transcripts,coverages,strands) if s == '+']
                            summed_coverage = sum([c for t,c in stranded_subset])
                            for t,c in stranded_subset:
                                p = c/summed_coverage
                                reads_to_assign = weight_triple[0]*p
                                assigned_reads[t] += reads_to_assign
                                for e in bridge.split('|'):
                                    if e in transcript_edges[t]:
                                        transcript_edges[t][e] += reads_to_assign
                    
                    if weight_triple[1]:
                        if strand_number.get('-',0) == 0:
                            weight_triple[1] = 0
                        elif strand_number.get('-',0) == 1:
                            assigned_reads[compatible_transcripts[which(strands,'-')[0]]] += weight_triple[1]
                        else:
                            stranded_subset = [(t,c) for t,c,s in zip(compatible_transcripts,coverages,strands) if s == '-']
                            summed_coverage = sum([c for t,c in stranded_subset])
                            for t,c in stranded_subset:
                                p = c/summed_coverage
                                reads_to_assign = weight_triple[1]*p
                                assigned_reads[t] += reads_to_assign
                                for e in bridge.split('|'):
                                    if e in transcript_edges[t]:
                                        transcript_edges[t][e] += reads_to_assign
                    
                    if weight_triple[2]:
                        summed_coverage = sum(coverages)
                        for t,c in zip(compatible_transcripts,coverages):
                            p = c/summed_coverage
                            reads_to_assign = weight_triple[2]*p
                            assigned_reads[t] += reads_to_assign
                            for e in bridge.split('|'):
                                if e in transcript_edges[t]:
                                    transcript_edges[t][e] += reads_to_assign
                    
                    unassigned_bridges[bridge] = [a-b for a,b in zip(unassigned_bridges[bridge],weight_triple)]
                    
                    
            # After all reads are assigned, update priors
            priors = [round(assigned_reads[t],digits) for t in transcripts]
            priors_dict = dict(zip(transcripts,priors))
            if any([i == 0 for i in priors]):
                delta = 0
                continue
            
            delta = sum([abs(a-b)/b for a,b in zip(priors,previous_priors)])/len(priors)
            # print(priors,delta)
            
        total_unassigned = sum([sum(i) for i in unassigned_bridges.values()])
        residuals = (total_unassigned/total_reads)+.01
        cumulative_deviation = 0
        deviation_list = []
        for t,edges in transcript_edges.items():
            coverage = priors_dict[t]/self.bridge_length(t)
            # Calculate deviance in coverage of all edges in the transcript
            grounded_edges = [v/self.bridge_length(k) for k,v in edges.items() if self.bridge_length(k)]
            deviation = sum([
                1/(1-(abs(edge_cov-coverage+0.001)/abs(edge_cov+coverage+0.001))+0.001)
                for edge_cov in grounded_edges
            ])/len(grounded_edges)
            deviation_list.append(deviation)
            cumulative_deviation += deviation
        
        # print(transcript_edges)
        locus_score = cumulative_deviation*residuals / len(transcripts) 
        output_dict = {}
            
        return locus_score, priors, deviation_list
    
    def optimal_transcripts(self,terminal_delta=1,maximum_in_group=8):
        """Performs a best-of-fit EM algorithm
        to iteratively improve the 'reads_assigned' values
        in the 'transcripts' dict. Returns a new dict
        with the optimal subset of transcripts and their
        read coverages."""
        optimal_transcripts = {}
        if not self.transcripts:
            self.maximum_flow()
        
        # Get overlap groups from the set of transcripts
        transcript_list = list(self.transcripts.keys())
        overlapping = overlap_groups(transcript_list)
        for o in overlapping:
            # If the number of models in o exceeds
            # maximum_in_group, remove the models with the
            # fewest reads assigned by maximum flow.
            if len(o) > maximum_in_group:
                o = [
                    o[j]
                    for j in order([
                        self.transcripts[i]['reads_assigned']
                        for i in o
                    ],reverse=True)[:maximum_in_group]
                ]
            
            best_score = None
            best_combo = []
            l = min([int(i.split(',')[0]) for i in o])
            r = max([int(i.split(',')[-2]) for i in o])
            transcript_combinations = combinations(o)
            for combo in transcript_combinations:
                # Calculate of Locus Score for each proposed
                # combo in the overlap group
                if not combo:
                    continue
                
                flow_scores = [self.transcripts[i]['reads_assigned'] for i in combo]
                locus_score, updated_scores, deviations = self.EM_read_assignment(l,r,combo,flow_scores,terminal_delta=terminal_delta)
                if best_score is None:
                    best_score = locus_score
                    best_combo = dict([(combo,{'reads_assigned':score,'deviation':dev,'locus_score':locus_score}) for combo,score,dev in zip(combo,updated_scores,deviations)])
                elif locus_score < best_score:
                    best_score = locus_score
                    best_combo = dict([(combo,{'reads_assigned':score,'deviation':dev,'locus_score':locus_score}) for combo,score,dev in zip(combo,updated_scores,deviations)])
            
            # After iterating over all combinations, store
            # add the highest scoring combo to optimal_transcripts
            optimal_transcripts.update(best_combo)
        
        # Make a modified copy of self.transcripts, only including
        # optimal_transcripts and their EM assigned reads
        for k in optimal_transcripts.keys():
            optimal_transcripts[k]['strand'] = self.transcripts[k]['strand']
        
        return optimal_transcripts
    
    def simplify(self):
        """Removes all trivial nodes from the network,
        i.e. all elements of the network that cannot be part of
        a transcript model. Returns a reduced copy of the object.
        """
        usable_nodes = []
        for i in range(self.number_of_nodes):
            if len(self.sources_from_sink(i)) > 0 and len(self.sinks_from_source(i)) > 0:
                if self.node_types[i] in ['D>','D<','A>','A<']:
                    if 'J' not in flatten(self.explode_edges(i)):
                        continue
                
                usable_nodes.append(i)
        
        # Find adjacent end features in usable_nodes and collapse them
        # into the most abundant signal
        last_node_type = ''
        last_score = 0
        redundant_ends = []
        for i in range(len(usable_nodes)):
            current_type = self.node_types[usable_nodes[i]]
            current_score = self.node_scores[usable_nodes[i]]
            if last_node_type == current_type and current_type in ['S>','S<','E>','E<']:
                if current_score > last_score:
                    redundant_ends.append(usable_nodes[i-1])
                elif current_score == last_score:
                    if current_type in ['S<','E>']:
                        redundant_ends.append(usable_nodes[i])
                    elif current_type in ['S>','E<']:
                        redundant_ends.append(usable_nodes[i-1])
                else:
                    redundant_ends.append(usable_nodes[i])
            
            last_node_type = current_type
            last_score = current_score
        
        usable_nodes = [i for i in usable_nodes if i not in redundant_ends]
        new_scores = [self.node_scores[i] for i in usable_nodes]
        
        new_daisychain = DaisyChain(
            chrom=self.chrom,
            start=self.start,
            node_locs=[self.node_locs[i] for i in usable_nodes],
            node_types=[self.node_types[i] for i in usable_nodes],
            node_scores=new_scores
        )
        return new_daisychain


