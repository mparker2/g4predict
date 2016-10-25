import sys
import os
import unittest
import regex

sys.path.append(
    os.path.abspath(os.path.dirname(
            os.path.dirname(
                __file__))))

import g4funcs as g4


class TestG4Regex(object):
    '''
    Some general tests for the G4Regex class.
    '''

    def test_tetrad_params(self):
        if self.test_params.get('tetrad_kwargs', False):
            self.assertEqual(self.g4regex._params['tetrad_kwargs']['start'],
                             self.test_params['tetrad_kwargs']['start'])
            self.assertEqual(self.g4regex._params['tetrad_kwargs']['stop'],
                             self.test_params['tetrad_kwargs']['stop'])

    def test_loop_params(self):
        if self.test_params.get('loop_kwargs_list', False):
            loops = [
                (x['start'],
                 x['stop']) for x in self.test_params['loop_kwargs_list']]
            for i, (start, stop) in enumerate(loops):
                self.assertEqual(
                    self.g4regex._params['loop_kwargs_list'][i]['start'], start
                )
                self.assertEqual(
                    self.g4regex._params['loop_kwargs_list'][i]['stop'], stop
                )

    def test_bulge_params(self):
        if self.test_params.get('bulge_kwargs', False):
            self.assertEqual(self.g4regex._params['bulge_kwargs']['start'],
                             self.test_params['bulge_kwargs']['start'])
            self.assertEqual(self.g4regex._params['bulge_kwargs']['stop'],
                             self.test_params['bulge_kwargs']['stop'])
            self.assertEqual(
                self.g4regex._params['bulge_kwargs']['bulges_allowed'],
                self.test_params['bulge_kwargs']['bulges_allowed'])

    def test_regex_production(self):
        self.assertIn(self.positive_strand_regex_example,
                      self.g4regex._regex['+'])
        self.assertIn(self.negative_strand_regex_example,
                      self.g4regex._regex['-'])

        self.assertEqual(len(self.g4regex._regex['+']),
                         self.len_regex)
        self.assertEqual(len(self.g4regex._regex['+']),
                         len(self.g4regex._regex['-']))

        if self.test_params.get('soft_mask', False):
            self.assertListEqual(self.g4regex._regex_flags, [])
        else:
            self.assertListEqual(self.g4regex._regex_flags,
                                 [regex.IGNORECASE])

    def test_matching(self):
        for seq, match in self.patterns:
            bed_12_gen = self.g4regex.get_g4s_as_bed(
                seq,
                seq_id='test',
                use_bed12=True)
            for m in match:
                self.assertEqual(next(bed_12_gen), m)
            with self.assertRaises(StopIteration):
                next(bed_12_gen)


class TestG4RegexTwoTetrad(TestG4Regex, unittest.TestCase):
    '''
    changing the default tetrad number only
    '''

    def setUp(self):
        self.test_params = dict(
            tetrad_kwargs=dict(start=2, stop=2),
        )
        self.g4regex = g4.G4Regex(**self.test_params)
        self.len_regex = 1
        self.positive_strand_regex_example = (
            '(?P<tet0>GG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GG)'
        )
        self.negative_strand_regex_example = (
            '(?P<tet0>CC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CC)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>CC)'
        )

        self.patterns = [
            # first pattern should have one match
            ['AAGGACTGGATGGTTTGGTTT',
             ['test\t2\t18\tPG4_2t0b\t28.0\t+\t2\t'
              '18\t85,118,209\t4\t2,2,2,2\t0,5,9,14']],
            # loops are too long here, should not match anything
            ['AAGGACTAAAAAATGGATGGTTTGGTTT',
             []],
            # soft mask is off, should still match
            ['AAGGACTggatggtttggTTT',
             ['test\t2\t18\tPG4_2t0b\t28.0\t+\t2\t'
              '18\t85,118,209\t4\t2,2,2,2\t0,5,9,14']]
        ]


class TestG4RegexTwoOrThreeTetrad(TestG4Regex, unittest.TestCase):
    '''
    allows two or three tetrad G Quadruplexes
    '''

    def setUp(self):
        self.test_params = dict(
            tetrad_kwargs=dict(start=2, stop=3),
        )
        self.g4regex = g4.G4Regex(**self.test_params)
        self.len_regex = 2
        # test three tetrad positive example...
        self.positive_strand_regex_example = (
            '(?P<tet0>GGG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GGG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GGG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GGG)'
        )
        # ...and two tetrad negative example
        self.negative_strand_regex_example = (
            '(?P<tet0>CC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CC)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>CC)'
        )

        self.patterns = [
            # first pattern should have one match
            ['AAGGACTGGATGGTTTGGTTT',
             ['test\t2\t18\tPG4_2t0b\t28.0\t+\t2\t'
              '18\t85,118,209\t4\t2,2,2,2\t0,5,9,14']],
            # matches one 3 tetrad and two 2 tetrad quadruplexes
            ['AAGGGACTGGGATGGGTTTGGGTTT',
             ['test\t2\t21\tPG4_2t0b\t23.5\t+\t2\t'
              '21\t85,118,209\t4\t2,2,2,2\t0,6,11,17',

              'test\t3\t21\tPG4_2t0b\t25.0\t+\t3\t'
              '21\t85,118,209\t4\t2,2,2,2\t0,5,10,16',

              'test\t2\t22\tPG4_3t0b\t48.0\t+\t2\t'
              '22\t85,118,209\t4\t3,3,3,3\t0,6,11,17']]
        ]


class TestG4RegexUnequalLoopLengths(TestG4Regex, unittest.TestCase):

    def setUp(self):
        self.test_params = dict(
            loop_kwargs_list=[
                dict(start=1, stop=x, allow_G=True) for x in [12, 7, 7]
            ],
        )
        self.g4regex = g4.G4Regex(**self.test_params)
        self.len_regex = 1
        self.positive_strand_regex_example = (
            '(?P<tet0>GGG)(?P<loop0>[ACGT]{1,12}?)'
            '(?P<tet1>GGG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GGG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GGG)'
        )
        self.negative_strand_regex_example = (
            '(?P<tet0>CCC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CCC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CCC)(?P<loop2>[ACGT]{1,12}?)'
            '(?P<tet3>CCC)'
        )

        self.patterns = [
            # should have one match
            ['AAGGGACTAAAAAATGGGATGGGTTTGGGTTT',
             ['test\t2\t29\tPG4_3t0b\t37.5\t+\t2\t'
              '29\t85,118,209\t4\t3,3,3,3\t0,13,18,24']],
            # if loop 3 is long, should not match anything,
            ['AAGGGACTTGGGATGGGTTAAAAAATGGGTTT',
             []],
            # should match the same as pattern 1 but on the negative strand
            ['AACCCACTTCCCATCCCTTAAAAAATCCCTTT',
             ['test\t2\t29\tPG4_3t0b\t37.5\t-\t2\t'
              '29\t85,118,209\t4\t3,3,3,3\t0,7,12,24']]
        ]


class TestG4RegexOneBulge(TestG4Regex, unittest.TestCase):

    def setUp(self):
        self.test_params = dict(
            bulge_kwargs=dict(bulges_allowed=1, start=1, stop=5)
        )
        self.g4regex = g4.G4Regex(**self.test_params)

        self.len_regex = 9
        self.positive_strand_regex_example = (
            '(?P<tet0>GGG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GGG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GGG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<btet3_1>G{2})'
            '(?P<bul3>[AT]{1,5})'
            '(?P<btet3_2>G{1})'
        )
        self.negative_strand_regex_example = (
            '(?P<tet0>CCC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CCC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CCC)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<btet3_1>C{2})'
            '(?P<bul3>[AT]{1,5})'
            '(?P<btet3_2>C{1})'
        )

        self.patterns = [
            # first pattern should have one match
            ['AAGGAGACTTGGGATGGGTTTGGGTTT',
             ['test\t2\t24\tPG4_3t1b\t40.0\t+\t'
              '2\t24\t85,118,209\t5\t2,1,3,3,3\t0,3,8,13,19']],
            # bulges with Cs in them are not allowed
            ['AAGGCAGACTTGGGATGGGTTTGGGTTT',
             []],
        ]


class TestG4RegexSoftMask(TestG4Regex, unittest.TestCase):

    def setUp(self):
        self.test_params = dict(
            soft_mask=True
        )
        self.g4regex = g4.G4Regex(**self.test_params)
        self.len_regex = 1
        self.positive_strand_regex_example = (
            '(?P<tet0>GGG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GGG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GGG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GGG)'
        )
        self.negative_strand_regex_example = (
            '(?P<tet0>CCC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CCC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CCC)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>CCC)'
        )

        self.patterns = [
            # first pattern should have one match
            ['AAGGGACTGGGATGGGTTTGGGTTT',
             ['test\t2\t22\tPG4_3t0b\t48.0\t+\t2\t'
              '22\t85,118,209\t4\t3,3,3,3\t0,6,11,17']],
            ['AAGGGACTgggatgggtttgggTTT',
             []],
        ]
