'''
Created on 27/06/2011

@author: dani
'''

from __future__ import division
from franklin.statistics import CachedArray, create_distribution
from os.path import join

from franklin.snv.snv_annotation import calculate_maf_frequency
from franklin.seq.readers import seqs_in_file

def create_pic_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None):
    'It creates the distribution of the pic'
    pics = CachedArray('f')
    for seq in seqs:
        for snv in seq.features:
            if 'pic' in snv.qualifiers:
                pics.append(snv.qualifiers['pic'])
    if list(pics):
        create_distribution(pics, labels={'title':'pic'},
                            distrib_fhand=distrib_fhand, bins=None,
                            plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

def create_het_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None):
    'It creates the distribution of the heterozygosity'
    hets = CachedArray('f')
    for seq in seqs:
        for snv in seq.features:
            if 'heterozygosity' in snv.qualifiers:
                hets.append(snv.qualifiers['heterozygosity'])
    if list(hets):
        create_distribution(hets, labels={'title':'heterozygosity'},
                            distrib_fhand=distrib_fhand, bins=None,
                            plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

def create_maf_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None):
    'It creates the distribution of the maf'
    mafs = CachedArray('f')
    for seq in seqs:
        for snv in seq.features:
            mafs.append(calculate_maf_frequency(snv))
    if list(mafs):
        create_distribution(mafs, labels={'title':'maf'},
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
