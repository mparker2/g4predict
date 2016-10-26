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


class TestClusterOverlapping(unittest.TestCase):

    def setUp(self):
        self.two_clusters_same_strand = StringIO(
            '1\t0\t50\ttest\t0\t+\n'
            '1\t25\t75\ttest\t0\t+\n'
            '1\t150\t200\ttest\t0\t+\n'
            '1\t175\t225\ttest\t0\t+\n'
        )
        self.two_clusters_opposite_strand = StringIO(
            '1\t0\t50\ttest\t0\t+\n'
            '1\t10\t60\ttest\t0\t+\n'
            '1\t15\t65\ttest\t0\t-\n'
            '1\t60\t110\ttest\t0\t-\n'
        )
        self.two_clusters_diff_chromosome = StringIO(
            '1\t0\t50\ttest\t0\t+\n'
            '1\t25\t75\ttest\t0\t+\n'
            '2\t0\t50\ttest\t0\t+\n'
            '2\t40\t90\ttest\t0\t+\n'
        )

    def test_cluster_overlapping(self):
        n_clusters_same_strand = 0
        for clust in g4.cluster_overlapping(self.two_clusters_same_strand):
            self.assertEqual(len(clust), 2)
            n_clusters_same_strand += 1
        self.assertEqual(n_clusters_same_strand, 2)

        n_clusters_opposite_strand = 0
        for clust in g4.cluster_overlapping(self.two_clusters_opposite_strand):
            self.assertEqual(len(clust), 2)
            n_clusters_opposite_strand += 1
        self.assertEqual(n_clusters_opposite_strand, 2)

        two_clusters_diff_chromosome = 0
        for clust in g4.cluster_overlapping(self.two_clusters_diff_chromosome):
            self.assertEqual(len(clust), 2)
            two_clusters_diff_chromosome += 1
        self.assertEqual(two_clusters_diff_chromosome, 2)


class TestOverlapMethod(object):
    '''
    general tests for g4funcs.filter_overlapping and g4funcs.merge_overlapping
    '''

    one_record_cluster = [
        ['1', 1, 100, 'test', 40, '+']
    ]
    two_record_cluster = [
        ['1', 0, 100, 'test', 40, '+'],
        ['1', 50, 150, 'test', 30, '+']
    ]
    three_record_cluster = [
        ['1', 0, 100, 'test', 40, '+'],
        ['1', 50, 150, 'test', 30, '+'],
        ['1', 120, 160, 'test', 20, '+']
    ]
    five_record_cluster = [
        ['1', 0, 100, 'test', 40, '+'],
        ['1', 50, 150, 'test', 30, '+'],
        ['1', 120, 160, 'test', 20, '+'],
        ['1', 155, 200, 'test', 20, '+'],
        ['1', 165, 300, 'test', 50, '+']
    ]

    def test_filter_method(self):

        def sort_bed_key(record):
            return int(record.split()[1])

        self.assertListEqual(sorted(self.method(self.one_record_cluster),
                                    key=sort_bed_key),
                             self.one_record_cluster_output)
        self.assertListEqual(sorted(self.method(self.two_record_cluster),
                                    key=sort_bed_key),
                             self.two_record_cluster_output)
        self.assertListEqual(sorted(self.method(self.three_record_cluster),
                                    key=sort_bed_key),
                             self.three_record_cluster_output)
        self.assertListEqual(sorted(self.method(self.five_record_cluster),
                                    key=sort_bed_key),
                             self.five_record_cluster_output)


class TestFilterOverlapping(TestOverlapMethod, unittest.TestCase):

    def setUp(self):
        self.method = g4.filter_overlapping
        self.one_record_cluster_output = ['1\t1\t100\ttest\t40\t+']
        self.two_record_cluster_output = ['1\t0\t100\ttest\t40\t+']
        self.three_record_cluster_output = ['1\t0\t100\ttest\t40\t+',
                                            '1\t120\t160\ttest\t20\t+']
        self.five_record_cluster_output = ['1\t0\t100\ttest\t40\t+',
                                           '1\t120\t160\ttest\t20\t+',
                                           '1\t165\t300\ttest\t50\t+']


class TestMergeOverlapping(TestOverlapMethod, unittest.TestCase):

    def setUp(self):
        self.method = g4.merge_overlapping
        self.one_record_cluster_output = ['1\t1\t100\tPG4_cluster\t1\t+']
        self.two_record_cluster_output = ['1\t0\t150\tPG4_cluster\t2\t+']
        self.three_record_cluster_output = ['1\t0\t160\tPG4_cluster\t3\t+']
        self.five_record_cluster_output = ['1\t0\t300\tPG4_cluster\t5\t+']


if __name__ == '__main__':
    unittest.main()
