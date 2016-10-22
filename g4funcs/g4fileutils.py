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

class FileWriter(object):

    def close(self):
        self.file.close()

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

    def write(self, bed_record):
        self.file.write('{}\n'.format(bed_record)

class BamWriter(FileWriter):

    def __init__(self, fn, chrom_sizes_file_name, is_faidx=True):

        if pysam is None:
            raise ImportError('pysam must be installed to write bam output')

        headers = self._get_headers(chrom_sizes_file_name, is_faidx)
        self.file = pysam.AlignmentFile(bam_file_name, 'wb', header=headers)
        self.fn = fn

    def _get_headers(self, chrom_sizes_file, is_faidx=True):
        headers = {
            'HD': {'VN': 1.0},
            'SQ': []
        }

        def parse_chromsizes_file(record):
            name, length = record.split()
            return name, length

        def parse_faidx_file(record):
            name, length, *_ = record.split()
            return name, length

        if is_faidx:
            parse_method = parse_faidx_file
        else:
            parse_method = parse_chromsizes_file

        with open(chrom_sizes_file) as f:
            for record in f:
                name, length = parse_method(record)
                headers['SQ'].append({'LN': int(length), 'SN': name})
        return headers   

    def write(self, bed_record):
        chrom, start, end, name, _, strand = bed_record.split()[:6]
        bam_record = pysam.AlignedSegment()
        bam_record.reference_name = chrom
        bam_record.query_name = name
        bam_record.reference_start = start
        bam_record.reference_length = end - start
        bam_record.mapping_quality = 255
        self.file.write(bam_record)

def sort_bed_file(unsorted_fn, sorted_fn=None):
    '''
    make temp or output bed file and write the sorted results to it.
    returns the name of the file containing sorted results.
    '''
    sorted_bed, sorted_fn = BedWriter(sorted_fn)
    subprocess.call(['sort', '-k1,1', '-k2,2n', unsorted_fn],
                    stdout=sorted_bed.file)
    sorted_bed.close()
    return sorted_fn

def choose_bam_or_bed(general_params):
    if general_params['write_bam']:
        faidx_file = '{}.fai'.format(general_params['fasta']
        if os.path.exists(faidx_file):
            output_file = g4.BamWriter(general_params['bed'], 
                                       faidx_file, 
                                       is_faidx=True)
        elif os.path.exists(general_params['contigs']):
            output_file = g4.BamWriter(general_params['bed'], 
                                       general_params['contigs'], 
                                       is_faidx=False)
        else:
            raise IOError('If writing bam, contigs file must be supplied,'
                          ' or fasta must be indexed with faidx')
    else:
        output_file = g4.BedWriter(general_params['bed'])
    return output_file

def apply_filter_method(file_handle, filter_method):
    temp_bed_file = get_bed_file()
    for cluster in cluster_overlapping(file_handle):
        for record in filter_method(cluster):
            temp_bed_file.write(record)
    temp_bed_file.close()
    sorted_file = sort_bed_file(temp_bed_file.fn)
    with open(sorted_file) as b:
        for line in b:
            yield b.strip()

