'''Some functions to do statistics on the sequences

Created on 10/07/2009
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import numpy
from itertools import tee
from biolib.biolib_utils import draw_scatter

PLOT_LABELS = {'masked_seq_distrib' :{
                                 'title': 'Masked sequence length distribution',
                                 'xlabel':'Masked length',
                                 'ylabel':'Number of sequences'},
                'seq_length_distrib':{'title':'Sequence length distribution',
                                      'xlabel':'Sequence length',
                                      'ylabel': 'Number of sequences'},
                'qual_distrib'      :{'title':'Sequence quality distribution',
                                      'xlabel':'quality',
                                      'ylabel':'Number of bp'}}


def _write_distribution(fhand, distribution, bin_edges):
    '''It writes the given distribution in a file.

    It requires a list with the distribution and a list with the bin_edges.
    len(bin_edges) == len(distribution) + 1
    '''
    for name, numbers in (('values', distribution), ('bin_edges', bin_edges)):
        fhand.write(name + ': ')
        fhand.write(','.join((str(item) for item in numbers)))
        fhand.write('\n')
    fhand.flush()

def _create_distribution(numbers, labels, distrib_fhand=None, plot_fhand=None,
                         range_=None):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    bins = 20
    #we do the distribution
    distrib, bin_edges = numpy.histogram(numbers, bins=bins, new=True,
                                         range=range_)
    #we write the output
    if distrib_fhand is not None:
        _write_distribution(distrib_fhand, distrib, bin_edges)
    #do we have to plot it?
    if plot_fhand is not None:
        draw_scatter(x_axe=bin_edges[:-1], y_axe=distrib,
                     title=labels['title'], xlabel=labels['xlabel'],
                     ylabel=labels['ylabel'],
                     fhand=plot_fhand)
    return {'distrib':distrib, 'bin_edges':bin_edges}


def _masked_seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                              range_=None):
    '''It calculates the masked length distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    lengths = []
    for seq in sequences:
        length = 0
        for letter in str(seq):
            if letter.islower():
                length += 1
        lengths.append(length)
    labels = PLOT_LABELS['masked_seq_distrib']
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand,
                                 range_=range_)

def _seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                       range_=None):
    '''It calculates the length distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    lengths = [len(seq) for seq in sequences]
    labels = PLOT_LABELS['seq_length_distrib']
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand,
                                range_=range_)

def _seq_qual_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                     range_=None):
    '''It calculates the qualities distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths

    qualities = []
    for seq in sequences:
        # This statistic can be run when there is no quality. In this cases it
        # should exit without warnings
        if seq.qual is None:
            return

        for value in seq.qual:
            qualities.append(value)
    labels = PLOT_LABELS['qual_distrib']
    return _create_distribution(qualities, labels, distrib_fhand, plot_fhand,
                                range_=range_)

#pylint:disable-msg=W0613
def seq_distrib(kind, sequences, distrib_fhand=None, plot_fhand=None,
                range_=None):
    'It returns the stat function depending on the stat it calculates(type)'
    stats = {'masked_seq_distrib'  :_masked_seq_length_distrib,
             'seq_length_distrib'  :_seq_length_distrib,
             'qual_distrib'        :_seq_qual_distrib }
    return stats[kind](sequences, distrib_fhand=distrib_fhand,
                       plot_fhand=plot_fhand, range_=range_)

def seq_distrib_diff(seqs1, seqs2, kind, distrib_fhand=None, plot_fhand=None):
    'It return the difference between different distributions of the given type'

    #get labesl depending on the stat type
    labels = PLOT_LABELS[kind]

    # I need two iterators for each ditribution calcule
    seqs1, seqs1_ = tee(seqs1, 2)
    seqs2, seqs2_ = tee(seqs2, 2)

    # first I need to calculate the distribution to calculate the combinen range
    distrib1 = seq_distrib(kind, seqs1)
    distrib2 = seq_distrib(kind, seqs2)

    range_ = _get_combined_distrib_range(distrib1, distrib2)

    #Now I calculate the distibution with the same range, in order to be able
    # to get the difference
    distrib1 = seq_distrib(kind, seqs1_, range_=range_)
    distrib2 = seq_distrib(kind, seqs2_, range_=range_)


    diff_distrib   = []
    diff_bin_edges = distrib1['bin_edges']
    # now a subtract distrib1 from distrib2
    for i in range(len(distrib1['distrib'])):
        diff = distrib1['distrib'][i] - distrib2['distrib'][i]
        diff_distrib.append(diff)

    #we write the output
    if distrib_fhand is not None:
        _write_distribution(distrib_fhand, diff_distrib, diff_bin_edges)
    #do we have to plot it?
    if plot_fhand is not None:
        draw_scatter(x_axe=diff_bin_edges[:-1], y_axe=diff_distrib,
                     title=labels['title'], xlabel=labels['xlabel'],
                     ylabel=labels['ylabel'],
                     fhand=plot_fhand)

    return {'distrib':diff_distrib, 'bin_edges':diff_bin_edges}



def _get_combined_distrib_range(distrib1, distrib2):
    'Giving two distribution return the range that intersecst both'

    min1 = min(distrib1['bin_edges'])
    max1 = max(distrib1['bin_edges'])
    min2 = min(distrib2['bin_edges'])
    max2 = max(distrib2['bin_edges'])

    if min1 > min2:
        min_ = min2
    else:
        min_ = min1
    if max1 > max2:
        max_ = max1
    else:
        max_ = max2
    return (min_, max_)

def general_seq_statistics(sequences, distrib_fhand=None, plot_fhand=None):
    'It counts the total length of the sum of all sequences'
    stats   = {}
    stats['seq_length']        = 0
    stats['num_sequences']       = 0
    stats['masked_length']  = 0
    seq_quality_average          = []

    first = True
    for seq in sequences:
        seq_len = len(seq)
        stats['num_sequences'] += 1
        stats['seq_length']  += seq_len

        # masked seq length calcule
        for nucleotide in seq.seq:
            if nucleotide.islower():
                stats['masked_length'] += 1
        # quality average calcule, I save seq length and quality average to
        # calculate more accurate toral average. need to check if seqs have
        # quality
        if first:
            first = False
            try:
                qual = seq.qual
                has_qual = True
            except AttributeError:
                has_qual = False

        if has_qual:
            qual = seq.qual
            qual_average = sum(qual)/float(seq_len)
            seq_quality_average.append((seq_len, qual_average))

    # length Average calcule
    stats['seq_length_average'] = stats['seq_length']/ \
                                                   float(stats['num_sequences'])

    if has_qual:
        total_qual = 0
        for seq_len, average in seq_quality_average:
            total_qual += seq_len * average
        stats['qual_average'] = total_qual / float(stats['seq_length'])


    # Write results to file if a file has been passed
    if distrib_fhand is not None:
        distrib_fhand.write('## Stats\n')
        for stat in stats:
            distrib_fhand.write('%s:%f\n' % (stat, float(stats[stat])))
    return stats




