'''
Created on 04/07/2011

@author: dani
'''
from franklin.seq.readers import seqs_in_file
from franklin.snv.snv_statistics import (create_pic_distribution,
                                         create_het_distribution,
                                         create_maf_distribution)
from os.path import join, exists
import os
from optparse import OptionParser
def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-s', '--inseqfile', dest='inseqfile',
                      help='input sequence file, (pickle)')

    parser.add_option('-d', '--dir_out', dest='dir_out',
                      help='Output directory')

    parser.add_option('-g', '--group_kind', dest='group_kind',
        help='Group kind to use (library(LB), platform(PL), sample(SM))')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.inseqfile:
        in_fpath = options.inseqfile
    else:
        parser.error('Input seq file mandatory')

    if options.dir_out:
        dir_out = options.dir_out
    else:
        dir_out = '.'

    if options.group_kind:
        group_kind = options.group_kind
        if group_kind not in ('LB', 'PL', 'SM'):
            parser.error('Group kind must be LB, PL or SM')
    else:

        group_kind = None

    return (in_fpath, dir_out, group_kind)


STAT_ANALYSIS = {'pic': create_pic_distribution, 'het': create_het_distribution,
                'maf': create_maf_distribution}


def main():
    'The main part'
    in_fpath, dir_out, group_kind = set_parameters()

    do_general_analysis(in_fpath, dir_out, group_kind)

def get_groups(fpath):
    groups = {'LB':[], 'PL':[], 'SM':[]}

    for seq in seqs_in_file(open(fpath)):
        for snv in seq.get_features(kind='snv'):
            read_groups = snv.qualifiers['read_groups']
            for tags in read_groups.values():
                for group_kind, group in tags.items():
                    if group not in groups[group_kind]:
                        groups[group_kind].append(group)
    return groups

def do_general_analysis(seq_path, dir_out, group_kind):
    groups = get_groups(seq_path)
    summary_fpath = join(dir_out, 'summary.all.txt')
    if exists(summary_fpath):
        os.remove(summary_fpath)
    sum_fhand = open(join(dir_out, 'summary.all.txt'), 'a')
    for analysis in STAT_ANALYSIS:
        seqs = seqs_in_file(open(seq_path, 'r'))
        dist_fhand = open(join(dir_out, analysis+'_distrib.all.dist'), 'w')
        svg_fhand = open(join(dir_out, analysis+'_distrib.all.svg'), 'w')
        STAT_ANALYSIS[analysis](seqs, distrib_fhand=dist_fhand,
                                plot_fhand=svg_fhand, summary_fhand=sum_fhand)
        dist_fhand.close()
        svg_fhand.close()

        if group_kind:
            for group in groups[group_kind]:
                seqs = seqs_in_file(open(seq_path, 'r'))
                dist_fhand = open(join(dir_out, analysis+'_distrib.%s.dist' % group) , 'w')
                svg_fhand = open(join(dir_out, analysis+'_distrib.%s.svg' % group), 'w')
                STAT_ANALYSIS[analysis](seqs, distrib_fhand=dist_fhand,
                                plot_fhand=svg_fhand, summary_fhand=sum_fhand,
                                group_kind=group_kind, groups=[group])
                dist_fhand.close()
                svg_fhand.close()
    sum_fhand.close()

if __name__ == '__main__':
    main()

