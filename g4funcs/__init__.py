import sys
from itertools import groupby
from .g4regex import *
from .g4filter import *

def parse_fasta(fasta):
    '''
    Parse a fasta file and yield the seq_id and sequence of each record
    '''
    def is_header(line):
        return line.startswith('>')
    
    if fasta == '-':
        f = sys.stdin
    else:
        f = open(fasta)
    
    fasta_it = groupby(f, is_header)
    for h, group in fasta_it:
        # take first word of fasta header as name, remove '>'
        header = next(group).split()[0][1:]
        _, seq_it = next(fasta_it)
        seq = ''.join(str(x).strip() for x in seq_it)
        yield (header, seq)
    
    f.close()