import os
import sys
import subprocess
from tempfile import mkstemp
from itertools import groupby
try:
    import pysam
except ImportError:
    pysam = None

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
        yield header, seq
    f.close()


class BedWriter(FileWriter):

    def __init__(self, fn=None):
        '''
        get a writable bed file object, can be temporary file, stdout, or 
        specific file
         '''
        if fn is None:
            fd, fn = mkstemp(suffix='.bed')
            self.file = os.fdopen(fd, 'w')
        elif fn == '-':
            self.file = sys.stdout
        else:
            self.file = open(fn, 'w')
        self.fn = fn

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.file.close()

    def close(self):
        self.file.close()

    def write(self, bed_record):
        self.file.write('{}\n'.format(bed_record)

def sort_bed_file(unsorted_fn, sorted_fn=None):
    '''
    make temp or output bed file and write the sorted results to it.
    returns the name of the file containing sorted results.
    '''
    with BedWriter(sorted_fn) as sorted_bed:
        subprocess.call(['sort', '-k1,1', '-k2,2n', unsorted_fn],
                         stdout=sorted_bed.file)
    return sorted_bed.fn

