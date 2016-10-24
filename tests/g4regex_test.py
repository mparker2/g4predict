import sys
import os
import unittest
import regex

sys.path.append(
    os.path.abspath(os.path.dirname(
            os.path.dirname(
                __file__))))

import g4funcs as g4


class TestG4Regex(unittest.TestCase):

    def setUp(self):
        self.g4regex = g4.G4Regex(**self.params)


class TestG4RegexTetrads(TestG4Regex):

    def __init__(self):
        super().__init__()
        self.params = dict(
            tetrad_kwargs=dict(start=2, stop=2),
            soft_mask=False
        )

    def test_correct_params(self):
        self.assertEqual(self.g4regex._params['tetrad_kwargs']['start'], 2)
        self.assertEqual(self.g4regex._params['tetrad_kwargs']['stop'], 2)
        self.assertEqual(self.g4regex._params['loop_kwargs'][0]['start'], 1)
        self.assertEqual(self.g4regex._params['loop_kwargs'][0]['stop'], 7)

    def test_regex(self):
        self.assertIn(
            '(?P<tet0>GG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GG)',
            self.g4regex._regex['+']
        )

        self.assertIn(
            '(?P<tet0>CC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CC)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>CC)',
            self.g4regex._regex['-']
        )

        self.assertEqual(len(self.g4regex._regex['+']), 1)
        self.assertEqual(self._regex_flags[0], regex.IGNORECASE)

    def test_matching(self):

        # test bed 6 output
        gen1 = self.g4regex.get_g4s_as_bed(
            'AAGGACTGGATGGTTTGGTTT',
            seq_id='test',
            use_bed12=False)
        self.assertEqual(
            next(gen1),
            'test\t2\t18\tPG4_2t0b\t28.0\t+')
        self.assertRaises(StopIteration, next(gen1))

        # test bed 12 output
        gen2 = self.g4regex.get_g4s_as_bed(
            'AAGGACTGGATGGTTTGGTTT',
            seq_id='test',
            use_bed12=True)
        self.assertEqual(
            next(gen2),
            'test\t2\t18\tPG4_2t0b\t28.0\t+'
            '\t2\t18\t85,118,209\t4\t2,2,2,2\t0,5,9,14')
        self.assertRaises(StopIteration, next(gen2))

        # loops are too long, should produce empty generator
        gen3 = self.g4regex.get_g4s_as_bed(
            'AAGGACTAAAAAATGGATGGTTTGGTTT',
            seq_id='test')
        self.assertRaises(StopIteration, next(gen3))


class TestG4RegexUnequalLoops(TestG4Regex):

    def __init__(self):
        super().__init()
        self.params = dict(
            loop_kwargs_list=[
                dict(start=1, stop=x, allow_G=True) for x in [12, 7, 7]
            ],
            soft_mask=True
        )

    def test_correct_params(self):
        for i, stop in enumerate([7, 12, 7]):
            self.assertEqual(
                self.g4regex._params['loop_kwargs'][i]['start'], 1)
            self.assertEqual(
                self.g4regex._params['loop_kwargs'][i]['stop'], stop)

    def test_regex(self):
        self.assertIn(
            '(?P<tet0>GGG)(?P<loop0>[ACGT]{1,12}?)'
            '(?P<tet1>GGG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GGG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<tet3>GGG)',
            self.g4regex._regex['+']
        )
        self.assertIn(
            '(?P<tet0>CCC)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>CCC)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>CCC)(?P<loop2>[ACGT]{1,12}?)'
            '(?P<tet3>CCC)',
            self.g4regex._regex['-']
        )
        self.assertEqual(len(self.g4regex._regex['+']), 1)
        self.assertEqual(self._regex_flags, [])

    def test_matching(self):

        # test bed 6 output
        gen1 = self.g4regex.get_g4s_as_bed(
            'AAGGGACTAAAAAATGGGATGGGTTTGGGTTT',
            seq_id='test',
            use_bed12=False)
        self.assertEqual(
            next(gen1),
            'test\t2\t18\tPG4_3t0b\t28.0\t+')
        self.assertRaises(StopIteration, next(gen1))

        # test bed 12 output
        gen2 = self.g4regex.get_g4s_as_bed(
            'AAGGACTAAAAAATGGATGGTTTGGGTTT',
            seq_id='test',
            use_bed12=True)
        self.assertEqual(
            next(gen2),
            'test\t2\t18\tPG4_3t0b\t28.0\t+'
            '\t2\t18\t85,118,209\t4\t2,2,2,2\t0,5,9,14')
        self.assertRaises(StopIteration, next(gen2))

        # loops are too long, should produce empty generator
        gen3 = self.g4regex.get_g4s_as_bed(
            'AAGGACTAAAAAATGGATGGTTTGGGTTT',
            seq_id='test')
        self.assertRaises(StopIteration, next(gen3))


class TestG4RegexBulges(TestG4Regex):

    def __init__(self):
        super().__init()
        self.params = dict(
            bulge_kwargs=dict(bulges_allowed=1, start=1, stop=5),
        )

    def test_correct_params(self):
        self.assertEqual(self.g4regex._params['bulge_kwargs']['start'], 1)
        self.assertEqual(self.g4regex._params['bulge_kwargs']['stop'], 5)
        self.assertEqual(
            self.g4regex._params['bulge_kwargs']['bulges_allowed'], 1)

    def test_regex(self):
        self.assertEqual(len(self.g4regex._regex['+']), 5)
        self.assertIn(
            '(?P<tet0>GG)(?P<loop0>[ACGT]{1,7}?)'
            '(?P<tet1>GG)(?P<loop1>[ACGT]{1,7}?)'
            '(?P<tet2>GG)(?P<loop2>[ACGT]{1,7}?)'
            '(?P<btet3_1>G{1})'
            '(?P<bul3>[AT]{1,5})'
            '(?P<btet3_2>G{1})',
            self.g4regex._regex['+'])

    def test_matching(self):
        gen1 = self.g4regex.get_g4s_as_bed(
            'AAGGAGACTTGGGATGGGTTTGGGTTT',
            seq_id='test',
            use_bed12=False)
        self.assertEqual(
            next(gen1),
            'test\t2\t24\tPG4_3t1b\t40.0\t+')
        self.assertRaises(StopIteration, next(gen1))

if __name__ == '__main__':
    unittest.main()
