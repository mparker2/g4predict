'''
Functions for handling fasta/bed files

author: Matthew Parker
'''


import os
import sys
import gzip
import signal
import subprocess
from tempfile import mkstemp
from itertools import groupby


class FileWrapper(object):

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.file.close()

    def close(self):
        self.file.close()


class FastaReader(FileWrapper):

    def __init__(self, fasta, decode_method=str):
        if isinstance(fasta, str):
            self._open_fasta(fasta)
        else:
            self.file = fasta
            self.decode_method = decode_method

    def parse_fasta(self):
        d = self.decode_method

        def is_header(line):
            return d(line).startswith('>')

        try:
            fasta_it = groupby(self.file, is_header)
        except TypeError:
            raise IOError('Object passed to FastaReader is not iterable')

        for h, group in fasta_it:
            # take first word of fasta header as name, remove '>'
            header = d(next(group)).split()[0][1:]
            _, seq_it = next(fasta_it)
            seq = ''.join(d(x).strip() for x in seq_it)
            yield header, seq

    def _open_fasta(self, fasta):
        if fasta == '-':
            self.file = sys.stdin
            self.decode_method = str
        elif os.path.splitext(fasta)[1] == '.gz':
            self.file = gzip.open(fasta)
            self.decode_method = bytes.decode
        else:
            self.file = open(fasta)
            self.decode_method = str


class BedWriter(FileWrapper):
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

    def write(self, bed_record):
        self.file.write('{}\n'.format(bed_record))


def sort_bed_file(unsorted_fn):
    '''
    sort a bed file using unix sort and yield sorted records in generator.
    '''
    def default_sigpipe():
        '''fixes some broken pipe behaviour'''
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    s = subprocess.Popen(['sort', '-k1,1', '-k2,2n', unsorted_fn],
                         stdout=subprocess.PIPE, preexec_fn=default_sigpipe)

    for record in s.stdout:
        yield record.decode().strip()
    s.stdout.close()
