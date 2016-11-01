'''
G4Regex: A set of classes for creating custom regular expressions to match
putative G Quadruplexes.

author: Matthew Parker
'''
from collections import defaultdict
from copy import copy, deepcopy
from itertools import product
from operator import itemgetter
import regex

# DEFAULT PARAMETERS:
# start and stop are inclusive
# no bulges allowed by default (bulged G4s are less thermostable/well
# characterised plus predicting them increases time/memory
# requirements greatly)

PARAMETERS = dict(
    tetrad_kwargs=dict(start=3, stop=3),
    loop_kwargs_list=[dict(start=1, stop=7, allow_G=True) for _ in range(3)],
    bulge_kwargs=dict(bulges_allowed=0, start=1, stop=5),
    score_kwargs=dict(tetrad_score_factor=20,
                      loop_pen_factor=1.5,
                      bulge_pen_factor=5),
    inter_kwargs=dict(start=2, stop=3),
    soft_mask=False
)

# REGEX BASES:
# regular expressions are built from these base strings using the parameters
# specified.
# In places, regular expression strings are formatted multiple times,
# or require curly brackets after formatting, hence the heavily escaped {{{}}}

# end result e.g. (?P<loop1>[ACGT]{1,7})
LOOP_BASE = (
    '(?P<loop{n}>'                # uses named loops (n increases 5'->3')
    '[ACGT]{{{start},{stop}}}?)'  # by default, loops can contain any base
)

# if specified loops, G4s with no Gs in loops can be matched
LOOP_BASE_NO_G = (
    '(?P<loop{n}>'
    # b is replaced by G or C depending on the strand
    '[AT{b}]{{{start},{stop}}})'
)

# end result e.g. (?P<tet1>GGG)
TETRAD_BASE = '(?P<tet{{n}}>{base})'  # uses named tetrads (n increases 5'->3')

BULGED_TETRAD_BASE = (
    # first half of tetrad, length of run is set by t1
    '(?P<btet{n}_1>{base}{{{t1}}})'
    # bulge, numbered the same as the tetrad, can contain only AT
    '(?P<bul{n}>[AT]{{{start},{stop}}})'
    # second half of tetrad, t1 + t2 == length of other non-bulged tetrads
    '(?P<btet{n}_2>{base}{{{t2}}})'
)


class G4Regex:

    '''
    Class for predicting G Quadruplexes.

    G Quadruplex tetrad length, loop length and bulge length parameters can
    be set to vary range of Quadruplexes which are predicted.
    '''

    def __init__(self, **kwargs):
        self._params = deepcopy(PARAMETERS)

        # update parameters
        if kwargs.get('tetrad_kwargs', False):
            self._params['tetrad_kwargs'].update(kwargs['tetrad_kwargs'])

        if kwargs.get('loop_kwargs_list', False):
            # shorten default parameters
            self._params['loop_kwargs_list'] = self._params[
                'loop_kwargs_list'][:len(kwargs['loop_kwargs_list'])]
            for i, new in enumerate(kwargs['loop_kwargs_list']):
                self._params['loop_kwargs_list'][i].update(new)

        if kwargs.get('bulge_kwargs', False):
            self._params['bulge_kwargs'].update(kwargs['bulge_kwargs'])

        if kwargs.get('score_kwargs', False):
            self._params['score_kwargs'].update(kwargs['score_kwargs'])

        # kwargs for PartialG4Regex which inherits from this class
        if kwargs.get('inter_kwargs', False):
            self._params['inter_kwargs'].update(kwargs['inter_kwargs'])

        if kwargs.get('soft_mask', False):
            self._params['soft_mask'] = True

        # use case insensitive matching if soft masking is turned off.
        if self._params['soft_mask']:
            self._regex_flags = []
        else:
            self._regex_flags = [regex.IGNORECASE]

        # self._regex stores generated regular expressions
        self._regex = defaultdict(list)

        self._build_g4_regex()

    def _build_g4_regex(self):

        # generate separate regex for each strand (G4s on neg strand are still
        # matched using the pos strand as a guide)
        for base, strand in (('G', '+'), ('C', '-')):

            # make copies of kwargs so that modifications do not affect other
            # strand
            loop_kwargs_c = deepcopy(self._params['loop_kwargs_list'])
            tetrad_kwargs_c = copy(self._params['tetrad_kwargs'])
            bulge_kwargs_c = copy(self._params['bulge_kwargs'])

            # generate loop regex
            loop_regex = []
            for i, kw in enumerate(loop_kwargs_c):

                # for each loop, check if G's are allowed.
                allow_G = kw.pop('allow_G', True)
                if allow_G:
                    loop_regex.append(
                        LOOP_BASE.format(n=i, **kw))
                else:
                    # We allow C in the loops, but not G
                    allowed_base = 'C' if strand == '+' else 'G'
                    loop_regex.append(
                        LOOP_BASE_NO_G.format(n=i, b=allowed_base, **kw))

            # reverse loops for opposite strand
            if strand == '-':
                loop_regex = loop_regex[::-1]

            # create individual regexes for each tetrad number
            # this is slower than backrefs when not allowing bulges,
            # but makes it much easier to extend the method to allow bulges
            for t in range(tetrad_kwargs_c['start'],
                           tetrad_kwargs_c['stop'] + 1):
                tet_regex = TETRAD_BASE.format(base=base * t)

                # maximum number of bulges that we can allow up to
                n_bulged = bulge_kwargs_c.pop('bulges_allowed', 0)

                # minimum number of unbulged tetrads can allow
                n_unbulged = 4 - n_bulged

                # the bulge can fall anywhere inside the tetrad, e.g.
                #    GATGGG
                #    GGATGG
                #    GGGATG
                # are all valid bulges for a 4 tetrad G4.
                # We use itertools product to create all possible combinations
                # of positions where tetrads could be bulged...
                # we use 0 to designate unbulged.
                bulge_combinations = product(range(0, t), repeat=4)

                # remove combinations which fall below min unbulged limit
                bulge_combinations = [
                    x for x in bulge_combinations if x.count(0) >= n_unbulged]

                # now we build a regex for each bulge combination
                for comb in bulge_combinations:
                    g4_regex = []
                    for i in range(4):  # iter over tetrads
                        if comb[i] == 0:

                            # if comb is 0 add unbulged tetrad, use i to name
                            g4_regex.append(tet_regex.format(n=i))
                        else:

                            # if comb != 0 we add a bulged tetrad
                            bulge_regex = BULGED_TETRAD_BASE.format(
                                base=base, n=i,
                                # we use the value of comb[i] for t1
                                t1=comb[i],
                                # and the total tetrad length minus it for t2
                                t2=t - comb[i],
                                **bulge_kwargs_c)
                            g4_regex.append(bulge_regex)

                        # append a loop after each tetrad
                        try:
                            g4_regex.append(loop_regex[i])
                        except IndexError:
                            # no loop after last tetrad
                            pass

                    # add to ever increasing dict of regexes
                    self._regex[strand].append(''.join(g4_regex))

    def get_g4s_as_bed(self, seq, seq_id='unknown', use_bed12=True):
        '''
        query a sequence for G4s using G4Regex. Pass a seq_id to get fully
        formatted bed records.
        Predicted loops/tetrad positional information can be retained using
        bed12 format.
        '''

        for strand in '+-':
            for r in self._regex[strand]:
                for m in regex.finditer(r,
                                        seq,
                                        overlapped=True,
                                        *self._regex_flags):
                    if use_bed12:
                        yield self._format_bed12(m, seq_id, strand)
                    else:
                        yield self._format_bed6(m, seq_id, strand)
                # clear re cache to save memory
                regex.purge()

    def _format_bed6(self, match, seq_id, strand):
        '''
        format a bed6 entry
        '''
        # use groupdict to count bulges and tetrads, to name the PG4
        gd = match.groupdict()
        tetrads = [v for k, v in gd.items() if k.startswith('tet')]
        l_tetrad = len(tetrads[0])  # length of each tetrad in bp

        loops = [gd['loop{}'.format(x)] for x in (0, 1, 2)]
        loops = ','.join(str(len(x)) for x in loops)

        bulges = [k for k, v in gd.items() if k.startswith('btet')]
        bulge_pos = set(k[4] for k in bulges)
        n_bulges = len(bulge_pos)
        bulge_flag = sum(2 ** int(f) for f in bulge_pos)

        start, end = match.span(0)
        name = '{}t{}b{}l'.format(l_tetrad, bulge_flag, loops)
        score = self._score_g4(l_tetrad, n_bulges, end - start)

        # format the bed record
        bed6 = '\t'.join(('{}',)*6)
        return bed6.format(seq_id, start, end, name, score, strand)

    def _format_bed12(self, match, seq_id, strand):
        '''
        format a bed12 entry
        '''
        # tetrads are always first and last matched groups with only one
        # other group between them: use [::2] to get their spans
        tetrad_spans = [
            match.span(x + 1) for x in range(len(match.groups()))][::2]
        start, end = match.span(0)

        # use groupdict to count bulges and tetrads, to name the PG4
        gd = match.groupdict()
        tetrads = [v for k, v in gd.items() if k.startswith('tet')]
        l_tetrad = len(tetrads[0])  # length of each tetrad in bp

        loops = [gd['loop{}'.format(x)] for x in (0, 1, 2)]
        loops = ','.join(str(len(x)) for x in loops)

        bulges = [k for k, v in gd.items() if k.startswith('btet')]
        bulge_pos = set(k[4] for k in bulges)
        n_bulges = len(bulge_pos)
        bulge_flag = sum(2 ** int(f) for f in bulge_pos)

        start, end = match.span(0)
        name = '{}t{}b{}l'.format(l_tetrad, bulge_flag, loops)
        score = self._score_g4(l_tetrad, n_bulges, end - start)
        rgb = '85,118,209'  # nice blue colour...

        # positional info
        # tetrads shown as blocks, loops+bulges as gaps
        block_count = 4 + n_bulges
        block_sizes = ','.join(str(y - x) for x, y in tetrad_spans)
        block_starts = ','.join(
            str(x - start) for x, _ in tetrad_spans)

        # format the bed record
        bed12 = '\t'.join(('{}',)*12)
        return bed12.format(
            seq_id, start, end, name, score, strand,
            start, end,  # thickStart/End the same as chromStart/End
            rgb, block_count, block_sizes, block_starts)

    def _score_g4(self, l_tetrad, n_bulge, length, n_tetrad=4):
        '''
        currently 'score' is just number of tetrads - total length of
        G4 loops and bulges, minus a gap penalty for bulges.
        '''
        score_params = self._params['score_kwargs']
        base_score = score_params['tetrad_score_factor'] * l_tetrad
        loop_pen = score_params['loop_pen_factor'] * (
            length - l_tetrad * n_tetrad)
        bulge_pen = score_params['bulge_pen_factor'] * n_bulge

        return base_score - loop_pen - bulge_pen


class PartialG4Regex(G4Regex):

    def _build_g4_regex(self):
        for base, strand in (('G', '+'), ('C', '-')):
            # make copies of kwargs so that modifications do not affect other
            # strand
            loop_kwargs_c = deepcopy(self._params['loop_kwargs_list'])
            tetrad_kwargs_c = copy(self._params['tetrad_kwargs'])
            inter_kwargs_c = copy(self._params['inter_kwargs'])

            # generate loop regex
            loop_regex = []
            for i, kw in enumerate(loop_kwargs_c):

                # for each loop, check if G's are allowed.
                allow_G = kw.pop('allow_G', True)
                if allow_G:
                    loop_regex.append(
                        LOOP_BASE.format(n=i, **kw))
                else:
                    # We allow C in the loops, but not G
                    allowed_base = 'C' if strand == '+' else 'G'
                    loop_regex.append(
                        LOOP_BASE_NO_G.format(n=i, b=allowed_base, **kw))

            # create individual regexes for each tetrad number
            for t in range(tetrad_kwargs_c['start'],
                           tetrad_kwargs_c['stop'] + 1):
                tet_regex = TETRAD_BASE.format(base=base * t)

                loop_regex_c = loop_regex[:t-1]
                # reverse loops for opposite strand
                if strand == '-':
                    loop_regex_c = loop_regex_c[::-1]

                g4_regex = ''
                # create regex for range of partial G4s.
                for i in range(inter_kwargs_c['stop']):
                    g4_regex += tet_regex.format(n=i)

                    if i in range(inter_kwargs_c['start'] - 1,
                                  inter_kwargs_c['stop']):
                        self._regex[strand].append(''.join(g4_regex))

                    # append a loop after each tetrad
                    try:
                        g4_regex += loop_regex_c[i]
                    except IndexError:
                        # no loop after last tetrad
                        break

    def _format_bed6(self, match, seq_id, strand):
        '''
        format a bed6 entry
        '''

        n_tetrad = match.re.pattern.count('tet')

        l_tetrad = len(match.group(1))  # length of each tetrad in bp

        start, end = match.span(0)
        name = 'PG4_{}t_{}'.format(l_tetrad, n_tetrad)
        score = self._score_g4(l_tetrad, 0, end - start, n_tetrad)

        # format the bed record
        bed6 = '\t'.join(('{}',)*6)
        return bed6.format(seq_id, start, end, name, score, strand)

    def _format_bed12(self, match, seq_id, strand):
        '''
        format a bed12 entry
        '''
        # tetrads are always first and last matched groups with only one
        # other group between them: use [::2] to get their spans
        tetrad_spans = [
            match.span(x + 1) for x in range(len(match.groups()))][::2]
        start, end = match.span(0)

        n_tetrad = match.re.pattern.count('tet')

        l_tetrad = len(match.group(1))  # length of each tetrad in bp

        name = 'PG4_{}t_{}'.format(l_tetrad, n_tetrad)
        score = self._score_g4(l_tetrad, 0, end - start, n_tetrad)
        rgb = '85,118,209'  # nice blue colour...

        # positional info
        # tetrads shown as blocks, loops+bulges as gaps
        block_count = n_tetrad
        block_sizes = ','.join(str(y - x) for x, y in tetrad_spans)
        block_starts = ','.join(
            str(x - start) for x, _ in tetrad_spans)

        # format the bed record
        bed12 = '\t'.join(('{}',)*12)
        return bed12.format(
            seq_id, start, end, name, score, strand,
            start, end,  # thickStart/End the same as chromStart/End
            rgb, block_count, block_sizes, block_starts)
