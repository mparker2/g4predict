import sys
import os
import unittest

sys.path.append(
    os.path.abspath(os.path.dirname(
            os.path.dirname(
                __file__))))

import g4funcs as g4


class TestG4OverlapMethod(object):
    '''
    test the filter overlapping method
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
        self.assertCountEqual(self.method(self.one_record_cluster),
                              self.one_record_cluster_output)
        self.assertCountEqual(self.method(self.two_record_cluster),
                              self.two_record_cluster_output)
        self.assertCountEqual(self.method(self.three_record_cluster),
                              self.three_record_cluster_output)
        self.assertCountEqual(self.method(self.five_record_cluster),
                              self.five_record_cluster_output)


class TestG4Filter(TestG4OverlapMethod, unittest.TestCase):

    def setUp(self):
        self.method = g4.filter_overlapping
        self.one_record_cluster_output = ['1\t1\t100\ttest\t40\t+']
        self.two_record_cluster_output = ['1\t0\t100\ttest\t40\t+']
        self.three_record_cluster_output = ['1\t0\t100\ttest\t40\t+',
                                            '1\t120\t160\ttest\t20\t+']
        self.five_record_cluster_output = ['1\t0\t100\ttest\t40\t+',
                                           '1\t120\t160\ttest\t20\t+',
                                           '1\t165\t300\ttest\t50\t+']


class TestG4Merge(TestG4OverlapMethod, unittest.TestCase):

    def setUp(self):
        self.method = g4.merge_overlapping
        self.one_record_cluster_output = ['1\t1\t100\tPG4_cluster\t1\t+']
        self.two_record_cluster_output = ['1\t0\t150\tPG4_cluster\t2\t+']
        self.three_record_cluster_output = ['1\t0\t160\tPG4_cluster\t3\t+']
        self.five_record_cluster_output = ['1\t0\t300\tPG4_cluster\t5\t+']


if __name__ == '__main__':
    unittest.main()
