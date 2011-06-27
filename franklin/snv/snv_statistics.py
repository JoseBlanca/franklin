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
            if 'pic' in snv.qualifiers['filters']:
                pics.append(snv.qualifiers['filters']['pic'][None])
    if list(pics):
        create_distribution(pics, labels=None, distrib_fhand=distrib_fhand,
                            bins=None, plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

def create_het_distribution(seqs, distrib_fhand=None, plot_fhand=None,
                            summary_fhand=None):
    'It creates the distribution of the heterozygosity'
    hets = CachedArray('f')
    for seq in seqs:
        for snv in seq.features:
            if 'heterozygosity' in snv.qualifiers['filters']:
                hets.append(snv.qualifiers['filters']['heterozygosity'][None])
    if list(hets):
        create_distribution(hets, labels=None, distrib_fhand=distrib_fhand,
                            bins=None, plot_fhand=plot_fhand, range_=None,
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
        create_distribution(mafs, labels=None, distrib_fhand=distrib_fhand,
                            bins=None, plot_fhand=plot_fhand, range_=None,
                            summary_fhand=summary_fhand, calculate_freqs=False,
                            remove_outliers=False)

#maf_distrib_fhand

STAT_ANALYSIS = {'pic': create_pic_distribution,
                'het': create_het_distribution,
                'maf': create_maf_distribution,
                }

def do_snv_stats(seq_path, out_dir):
    'It performs snv statistics'
    for analysis in STAT_ANALYSIS:
        seqs = seqs_in_file(open(seq_path.last_version, 'r'))
        dist_fhand = open(join(out_dir, analysis+'_distrib.dist'), 'w')
        svg_fhand = open(join(out_dir, analysis+'_distrib.svg'), 'w')
        sum_fhand = open(join(out_dir, analysis+'_distrib.sum'), 'w')
        STAT_ANALYSIS[analysis](seqs, distrib_fhand=dist_fhand,
                                plot_fhand=svg_fhand, summary_fhand=sum_fhand)








