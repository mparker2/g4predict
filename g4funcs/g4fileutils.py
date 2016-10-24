'''
Functions for handling fasta/bed files

author: Matthew Parker
'''


import os
import sys
import gzip
import subprocess
from tempfile import mkstemp
from itertools import groupby


def parse_fasta(fasta):
    '''
    Parse a fasta file and yield the seq_id and sequence of each record
    '''

    if fasta == '-':
        f = sys.stdin
        decode_method = str
    elif os.path.splitext(fasta)[1] == '.gz':
        f = gzip.open(fasta)
        decode_method = bytes.decode
    else:
        f = open(fasta)
        decode_method = str

    def is_header(line):
        return decode_method(line).startswith('>')

    fasta_it = groupby(f, is_header)
    for h, group in fasta_it:
        # take first word of fasta header as name, remove '>'
        header = decode_method(next(group)).split()[0][1:]
        _, seq_it = next(fasta_it)
        seq = ''.join(decode_method(x).strip() for x in seq_it)
        yield header, seq
    f.close()


class BedWriter(object):
    '''
    get a writable bed file object, can be temporary file, stdout, or
    specific file. Default is to make new temp file in /tmp/
     '''

    def __init__(self, fn=None):
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
        self.file.write('{}\n'.format(bed_record))


def sort_bed_file(unsorted_fn):
    '''
    sort a bed file using unix sort and yield sorted records in generator.
    '''
    s = subprocess.Popen(['sort', '-k1,1', '-k2,2n', unsorted_fn],
                         stdout=subprocess.PIPE)
    for record in s.stdout:
        yield record.decode().strip()
