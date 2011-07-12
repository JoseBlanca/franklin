'''
Created on 27/06/2011

@author: dani
'''

from __future__ import division
from franklin.statistics import CachedArray, create_distribution
from os.path import join

from franklin.snv.snv_annotation import (calculate_maf_frequency,
                                         calculate_heterozygosity,
                                         calculate_pic)
from franklin.seq.readers import seqs_in_file

def create_pic_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None, read_groups=None,
                            group_kind=None, groups=None):
    'It creates the distribution of the pic (not takes in account ref allele)'
    title = 'pic'
    if groups and group_kind:
        title = 'pic (%s: %s)' % (group_kind, ','.join(groups))

    pics = CachedArray('f')
    for seq in seqs:
        for snv in seq.get_features('snv'):
            if not group_kind and 'pic' in snv.qualifiers:
                pic = snv.qualifiers['pic']
            else:
                pic = calculate_pic(snv, group_kind=group_kind, groups=groups)
            if pic is not None:
                pics.append(pic)
    if list(pics):
        create_distribution(pics, labels={'title':title},
                            distrib_fhand=distrib_fhand, bins=None,
                            plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

def create_het_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None, group_kind=None,
                            groups=None, ploidy=2):
    '''It creates the distribution of the heterozygosity
    (not takes in account ref allele)'''
    title = 'heterozygosity'
    if groups and group_kind:
        title = 'heterozygosity (%s: %s)' % (group_kind, ','.join(groups))

    hets = CachedArray('f')
    for seq in seqs:
        for snv in seq.get_features('snv'):
            if not group_kind and 'heterozygosity' in snv.qualifiers:
                het = snv.qualifiers['heterozygosity']
            else:
                het = calculate_heterozygosity(snv, ploidy,
                                               group_kind=group_kind,
                                               groups=groups)
            if het is not None:
                hets.append(het)
    if list(hets):
        create_distribution(hets, labels={'title':title},
                            distrib_fhand=distrib_fhand, bins=None,
                            plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

def create_maf_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None, groups=None, group_kind=None):
    'It creates the distribution of the maf (not takes in account ref allele)'
    title = 'maf'
    if groups and group_kind:
        title = 'maf (%s: %s)' % (group_kind, ','.join(groups))

    mafs = CachedArray('f')
    for seq in seqs:
        for snv in seq.get_features('snv'):
            maf = calculate_maf_frequency(snv, groups=groups,
                                          group_kind=group_kind)
            if maf:
                mafs.append(maf)
    if list(mafs):
        create_distribution(mafs, labels={'title':title},
                            distrib_fhand=distrib_fhand, bins=None,
                            plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

STAT_ANALYSIS = {'pic': create_pic_distribution, 'het': create_het_distribution,
                'maf': create_maf_distribution}

def do_snv_stats(seq_path, out_dir):
    'It performs snv statistics'
    first_time = True
    for analysis in STAT_ANALYSIS:
        seqs = seqs_in_file(open(seq_path.last_version, 'r'))
        dist_fhand = open(join(out_dir, analysis+'_distrib.dist'), 'w')
        svg_fhand = open(join(out_dir, analysis+'_distrib.svg'), 'w')
        if first_time == True:
            sum_fhand = open(join(out_dir, 'distrib.sum'), 'w')
            first_time = False
        else:
            sum_fhand = open(join(out_dir, 'distrib.sum'), 'a')
        STAT_ANALYSIS[analysis](seqs, distrib_fhand=dist_fhand,
                                plot_fhand=svg_fhand, summary_fhand=sum_fhand)
        dist_fhand.close()
        svg_fhand.close()
        sum_fhand.close()
