'''
Created on 08/07/2011

@author: dani
'''

from optparse import OptionParser
from franklin.gff import GffFile
from itertools import ifilter
from operator import itemgetter
import sys

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--in_gff', dest='in_gff',
                      help='Input ggf file')

    parser.add_option('-o', '--out_bed', dest='out_bed',
                      help='Output bed file')

    parser.add_option('-k', '--feat_kinds', dest='feat_kinds',
                      help='Feature kind')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.in_gff is None:
        parser.error('Input seq file mandatory')
    else:
        in_fpath = options.in_gff

    if options.out_bed is None:
        out_fhand = sys.stdout
    else:
        out_fhand = open(options.out_bed, 'w')

    if options.feat_kinds:
        feat_kinds = options.feat_kinds
    else:
        parser.error('One feature kind is needed at least')

    return in_fpath, out_fhand, feat_kinds


def gff2bed(in_fpath, out_bed_fhand, feat_kinds=None, sort=True):
    'It converts gff file into bed file'
    gff = GffFile(in_fpath)
    feats = gff.features
    if feat_kinds:
        feats = ifilter(lambda x: x['type'] in feat_kinds,  feats)

    feats = ((feat['seqid'], feat['start'], feat['end']) for feat in feats)

    if sort:
        feats = sorted(feats, key=itemgetter(1))
        feats = sorted(feats, key=itemgetter(0))

    for feat in feats:
        feat = map(str, feat)
        line = '\t'.join(feat) + '\n'
        out_bed_fhand.write(line)
    out_bed_fhand.flush()

def main():
    'The main part'
    in_fpath, out_bed_fhand, feat_kinds = set_parameters()
    gff2bed(in_fpath, out_bed_fhand, feat_kinds)

if __name__ == '__main__':
    main()
