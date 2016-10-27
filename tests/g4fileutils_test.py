import sys
import os
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

sys.path.append(
    os.path.abspath(os.path.dirname(
            os.path.dirname(
                __file__))))

import g4funcs as g4


class TestParseFasta(unittest.TestCase):

    def setUp(self):
        self.fasta_file = g4.FastaReader(StringIO('''
>1
ATGTGTGAGTGTGTGAGTGTGTGTTTAGTGTTGAGAGT
AGTAGTGAGTGTGCGCGGGAGAGTGTCGGTAGTGTGTG
GATAGATAAGCTACGCATCAGCACATTATTATATATAT
GGGGCAGCGAGCAGCAGTCAGCATAGCATCAGCATCAG
AGTAGTCAG
>2
GGCATCGACTACGATGACGACATCGACTATCAGCATAC
GACTAGCATCAGCAGCTACGCTAGGACTAGATCAGCAT
CGTGCACGACTGATAGACTACGCTACGATCAGCT
>3test
GCATCAGCATCGATACTATTTATTATATATATTATATA
GATAGCTATCTACTATATCATTAATATATT
'''.strip()))
        self.seqs = [
            'ATGTGTGAGTGTGTGAGTGTGTGTTTAGTGTTGAGAGT'
            'AGTAGTGAGTGTGCGCGGGAGAGTGTCGGTAGTGTGTG'
            'GATAGATAAGCTACGCATCAGCACATTATTATATATAT'
            'GGGGCAGCGAGCAGCAGTCAGCATAGCATCAGCATCAG'
            'AGTAGTCAG',
            'GGCATCGACTACGATGACGACATCGACTATCAGCATAC'
            'GACTAGCATCAGCAGCTACGCTAGGACTAGATCAGCAT'
            'CGTGCACGACTGATAGACTACGCTACGATCAGCT',
            'GCATCAGCATCGATACTATTTATTATATATATTATATA'
            'GATAGCTATCTACTATATCATTAATATATT'
        ]
        self.seq_ids = ['1', '2', '3test']

    def test_parse_fasta(self):
        fasta_iter = self.fasta_file.parse_fasta()
        for test_id, test_seq, in zip(self.seq_ids, self.seqs):
            parsed_id, parsed_seq = next(fasta_iter)
            self.assertEqual(test_id, parsed_id)
            self.assertEqual(test_seq, parsed_seq)
        with self.assertRaises(StopIteration):
            next(fasta_iter)


class TestSortBed(unittest.TestCase):

    def setUp(self):
        unsorted_bed = [
            u'1\t100\t200\ttest\t10\t-',
            u'1\t0\t100\ttest\t10\t+',
            u'1\t150\t250\ttest\t10\t-',
            u'2\t150\t250\ttest\t10\t-',
            u'1\t120\t220\ttest\t10\t+',
            u'1\t10\t110\ttest\t10\t-',
            u'2\t120\t220\ttest\t10\t-',
        ]
        with g4.BedWriter() as bed:
            for r in unsorted_bed:
                bed.write(r)
        self.unsorted_bed_fn = bed.fn
        self.sorted_bed = [
            u'1\t0\t100\ttest\t10\t+',
            u'1\t10\t110\ttest\t10\t-',
            u'1\t100\t200\ttest\t10\t-',
            u'1\t120\t220\ttest\t10\t+',
            u'1\t150\t250\ttest\t10\t-',
            u'2\t120\t220\ttest\t10\t-',
            u'2\t150\t250\ttest\t10\t-',
        ]

    def test_sort_bed(self):
        sorted_output = list(g4.sort_bed_file(self.unsorted_bed_fn))
        self.assertEqual(sorted_output, self.sorted_bed)
