'''
Predict putative G Quadruplexes using an extension of the Quadparser method;
NB: Output from g4predict is unlikely to be correctly sorted. use unix sort.
Author: Matthew Parker;
'''

import os
import sys
from tempfile import mkstemp
import argparse
import subprocess
import g4funcs as g4



def parse_args():
    '''
    Get command line arguments
    '''
    
    def intra(args):
        '''
        parse inter args correctly and return G4Regex
        instance and general_params
        '''
        
        # deal with loop arguments
        for l_arg in ('min_loop', 'max_loop', 'allow_g'):
            v = args[l_arg].split(',')
            if len(v) == 1:
                v *= 3
            elif len(v) != 3:
                raise argparse.ArgumentTypeError(
                    '{} should be a single int or 3 comma'
                    ' separated ints'.format(l_arg))
            args[l_arg] = [int(x) for x in v]
        
        g4_params = dict(
            tetrad_kwargs = dict(
                start=args.pop('min_tetrad'),
                stop=args.pop('max_tetrad')),
            loop_kwargs_list = [
                dict(start=x, stop=y, allow_G=z) for x, y, z in zip(
                    args.pop('min_loop'),
                    args.pop('max_loop'),
                    args.pop('allow_g'))],
            bulge_kwargs = dict(
                bulges_allowed = args.pop('bulges'),
                start = args.pop('min_bulge'),
                stop = args.pop('max_bulge')),
            score_kwargs = dict(
                tetrad_score_factor = args.pop('tetrad_score_factor'),
                loop_pen_factor = args.pop('loop_pen_factor'),
                bulge_pen_factor = args.pop('bulge_pen_factor'))
            )
        
        return args, g4.G4Regex(**g4_params)
    
    def inter(args):
        '''
        parse intra args correctly and return PartialG4Regex
        instance and general_params
        '''
        
        # deal with loop arguments
        for l_arg in ('min_loop', 'max_loop', 'allow_g'):
            v = args[l_arg].split(',')
            if len(v) == 1:
                v *= args['max_g_runs']
            elif len(v) != (args['max_g_runs'] - 1):
                raise argparse.ArgumentTypeError(
                    '{} should be a single int or comma separated ints '
                    'of length --max-g-runs minus one'.format(l_arg))
            args[l_arg] = [int(x) for x in v]
        
        g4_params = dict(
            tetrad_kwargs = dict(
                start=args.pop('min_tetrad'),
                stop=args.pop('max_tetrad')),
            loop_kwargs_list = [
                dict(start=x, stop=y, allow_G=z) for x, y, z in zip(
                    args.pop('min_loop'),
                    args.pop('max_loop'),
                    args.pop('allow_g'))],
            inter_kwargs = dict(
                start = args.pop('min_g_runs'),
                stop = args.pop('max_g_runs')),
            score_kwargs = dict(
                tetrad_score_factor = args.pop('tetrad_score_factor'),
                loop_pen_factor = args.pop('loop_pen_factor'))
            )
        
        return args, g4.PartialG4Regex(**g4_params)    
    
    
    a = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # two subparsers, 'intra' and 'inter'
    sub = a.add_subparsers()
    intra_parser = sub.add_parser('intra', help='''
Predict complete, intramolecular PG4s (i.e. PG4s which form from one 
DNA/RNA strand). Uses the general pattern G{x}([ATGC]{y,z}G{x}){3}.
''')
    intra_parser.set_defaults(func=intra)    
    
    inter_parser = sub.add_parser('inter', help='''
Predict partial, intermolecular PG4s, which cannot form on there own but
require at least 1 other partial G4 from a different DNA/RNA molecule.
''')
    inter_parser.set_defaults(func=inter)    
    
    for p in (inter_parser, intra_parser):
        
        general = p.add_argument_group('General')
        general.add_argument(
            '-f', '--fasta', type=str, required=True,
            help='Input fasta file, use \'-\' to read from stdin')
        general.add_argument(
            '-b', '--bed', type=str, required=True,
            help='Output bed file, use \'-\' to write to stdout')
        general.add_argument(
            '-t', '--write-bed12', action='store_true', required=False,
            default=True,
            help='write bed12 output')
        general.add_argument(
            '-s', '--write-bed6', action='store_true', required=False,
            help='write bed6 output instead of bed12 (some information is lost)')
        general.add_argument(
            '-F', '--filter-overlapping', action='store_true',
            required=False, default=False,
            help='''
use filtering method to remove overlapping PG4s, yields the maximum number of
high scoring, non-overlapping PG4s
''')
        general.add_argument(
            '-M', '--merge-overlapping', action='store_true',
            required=False, default=False,
            help='''
use merge method to flatten overlapping PG4s, output is in bed6 and overrides
the --write-bed12 flag
''')

        score = p.add_argument_group('Score', description='''
Score method parameters. Scoring method is x*T - y*L - z*B, where T is the  
number of tetrads, L is the total length of all loops and bulges, B is the
number of bulges, and x, y, and z are parameters which can be set by the user. 
''')
        score.add_argument(
            '-x', '--tetrad-score-factor', type=float, required=False,
            default=20,
            help='factor to be multiplied by tetrad number in scoring')    
        score.add_argument(
            '-y', '--loop-pen-factor', type=float, required=False, default=1.5,
            help='factor to be multiplied by loop length in scoring')        

        tetrad = p.add_argument_group('Tetrads')
        tetrad.add_argument(
            '-tmin', '--min-tetrad', type=int, required=False, default=3,
            help='min number of tetrads allowed in predicted G4')
        tetrad.add_argument(
            '-tmax', '--max-tetrad', type=int, required=False, default=3,
            help='max number of tetrads in predicted G4')
    
        loop = p.add_argument_group('Loops')
        loop.add_argument(
            '-lmin', '--min-loop', type=str, required=False, default='1',
            help='''
min loop length of predicted G4, either a single int if all 
loops mins to the same length, or 3 comma separated ints for loops 1 up to 
3 (5'-> 3')
''')
        loop.add_argument(
            '-lmax', '--max-loop', type=str, required=False, default='7',
            help='''
max loop length of predicted G4, either a single int if all 
loop maxes to be the same, or 3 comma separated ints for loops 1 up to 
3 (5'-> 3')
''')
        loop.add_argument(
            '-G', '--allow-g', type=str, required=False, default='1',
            help='''
allow G in PG4 loops, use 0 to disallow G in all loops, or comma separated 0s 
or 1s to disallow G in specific loops
''')
    
    # Intra specific parser
    bulge = intra_parser.add_argument_group('Bulges (Intra Only)')
    bulge.add_argument(
        '-B', '--bulges', type=int, required=False, default=0,
        help='how many bulges to allow in PG4s, can have up to one per tetrad')    
    bulge.add_argument(
        '-bmin', '--min-bulge', type=int, required=False, default=1,
        help='min bulge length allowed in predicted G4')
    bulge.add_argument(
        '-bmax', '--max-bulge', type=int, required=False, default=5,
        help='min bulge length allowed in predicted G4')
    # this score parameter is only added to intra options
    score.add_argument(
        '-z', '--bulge-pen-factor', type=float, required=False, default=5,
        help='factor to be multiplied by bulge number in scoring')    
    
    # Inter specific parser
    inter = inter_parser.add_argument_group('Inter')
    inter.add_argument(
        '-rmin', '--min-g-runs', type=int, required=False, default=2,
        help='min runs of G to use to predict partial PG4s')
    inter.add_argument(
        '-rmax', '--max-g-runs', type=int, required=False, default=3,
        help='max runs of G to use to predict partial PG4s')   
    
    args = a.parse_args()
    if args.write_bed12 is args.write_bed6:
        a.error('--write-bed12 and --write-bed6 are mutually exclusive')
    elif args.write_bed6:
        args.write_bed12 = False
    
    if args.filter_overlapping is args.merge_overlapping is True:
        a.error(
            '--filter-overlapping and --merge-overlapping'
            ' are mutually exclusive')
    
    return args.func(vars(args))

if __name__ == '__main__':
    general_params, g4regex = parse_args()
    fasta = g4.parse_fasta(general_params['fasta'])
    
    # if we want to filter overlapping records, we need to write to file, then
    # sort the file before we do the filtering
    if general_params['filter_overlapping'] or general_params['merge_overlapping']:
        fd1, fn1 = mkstemp(suffix='_g4pred.bed')
        bed = os.fdopen(fd1, 'w')
    else:
        if general_params['bed'] != '-':
            bed = open(general_params['bed'], 'w')
        else:
            bed = sys.stdout
    
    
    for seq_id, seq in fasta:
        for record in g4regex.get_g4s_as_bed(seq, seq_id=seq_id,
                                             use_bed12=general_params['write_bed12']):
            bed.write('{}\n'.format(record))
    
    bed.close()
    
    if general_params['filter_overlapping'] or general_params['merge_overlapping']:
        
        # make temp bed file and write the sorted results to it
        fd2, fn2 = mkstemp(suffix='_g4sort.bed.')
        sorted_bed = os.fdopen(fd2)
        subprocess.call(['sort', '-k1,1', '-k2,2n', fn1],
                        stdout=sorted_bed)
        sorted_bed.close()
        
        # reopen sorted bedfile and filter it:
        if general_params['bed'] != '-':
            outbed = open(general_params['bed'], 'w')
        else:
            outbed = sys.stdout
        with open(fn2) as inbed:
            for cluster in g4.cluster_overlapping(inbed):
                
                # apply filter method...
                if general_params['filter_overlapping']:
                    records = g4.filter_overlapping(cluster)
                    for record in records:
                        record = '\t'.join([str(x) for x in record])
                        outbed.write('{}\n'.format(record))
                
                # ...or merge method
                elif general_params['merge_overlapping']:
                    record = g4.merge_overlapping(cluster)
                    record = '\t'.join([str(x) for x in record])
                    outbed.write('{}\n'.format(record))
        
        outbed.close()
        # clean up temp files
        os.remove(fn1)
        os.remove(fn2)