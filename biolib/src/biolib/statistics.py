'''Some functions to do statistics on sequences.'''

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
from biolib.biolib_utils import draw_scatter, FileCachedList


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

def _create_distribution(numbers, labels=None, distrib_fhand=None,
                         plot_fhand=None, range_=None, low_memory=True):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    bins = 20
    if labels is None:
        labels = {'title':'histogram', 'xlabel':'values', 'ylabel':'count'}
    #we do the distribution
    if low_memory:
        distrib, bin_edges = histogram(numbers, bins=bins, range_=range_)
    else:
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

def _masked_sequence_lengths(sequences, low_memory):
    'It returns a list with the sequences lengths'
    if low_memory:
        lengths = FileCachedList(float)
    else:
        lengths = []
    for seq in sequences:
        length = 0
        for letter in str(seq):
            if letter.islower():
                length += 1
        lengths.append(length)
    return lengths

def _sequence_lengths(sequences, low_memory):
    'It returns a list with the sequences lengths'
    if low_memory:
        lengths = FileCachedList(float)
    else:
        lengths = []
    for seq in sequences:
        lengths.append(len(seq))
    return lengths

def _sequence_qualitities(sequences, low_memory):
    'It returns a list with the sequence qualities'
    if low_memory:
        qualities = FileCachedList(float)
    else:
        qualities = []
    for seq in sequences:
        # This statistic can be run when there is no quality. In this cases it
        # should exit without warnings
        if seq.qual is None:
            continue
        for value in seq.qual:
            qualities.append(value)
    return qualities

def _masked_seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                              range_=None, low_memory=True):
    '''It calculates the masked length distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    lengths = _masked_sequence_lengths(sequences, low_memory)
    labels = PLOT_LABELS['masked_seq_distrib']
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand,
                                 range_=range_, low_memory=low_memory)

def _seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                       range_=None, low_memory=True):
    '''It calculates the length distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    lengths = _sequence_lengths(sequences, low_memory)
    labels = PLOT_LABELS['seq_length_distrib']
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand,
                                range_=range_, low_memory=low_memory)

def _seq_qual_distrib(sequences, distrib_fhand=None, plot_fhand=None,
                     range_=None, low_memory=True):
    '''It calculates the qualities distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    qualities = _sequence_qualitities(sequences, low_memory)
    labels = PLOT_LABELS['qual_distrib']
    return _create_distribution(qualities, labels, distrib_fhand, plot_fhand,
                                range_=range_, low_memory=low_memory)

#pylint:disable-msg=W0613
def seq_distrib(kind, sequences, distrib_fhand=None, plot_fhand=None,
                range_=None):
    'It returns the stat function depending on the stat it calculates(type)'
    stats = {'masked_seq_distrib'  :_masked_seq_length_distrib,
             'seq_length_distrib'  :_seq_length_distrib,
             'qual_distrib'        :_seq_qual_distrib }
    return stats[kind](sequences, distrib_fhand=distrib_fhand,
                       plot_fhand=plot_fhand, range_=range_)

def seq_distrib_diff(seqs1, seqs2, kind, distrib_fhand=None, plot_fhand=None,
                     low_memory=True):
    'It return the difference between different distributions of the given type'

    #get labesl depending on the stat type
    labels = PLOT_LABELS[kind]

    #the values for boths sequences
    values_functs = {
             'masked_seq_distrib': _masked_sequence_lengths,
             'seq_length_distrib': _sequence_lengths,
             'qual_distrib'      : _sequence_qualitities}
    values1 = values_functs[kind](seqs1, low_memory=low_memory)
    values2 = values_functs[kind](seqs2, low_memory=low_memory)

    #the range
    if low_memory:
        num_iters1, num_iters2 = values1.items(), values2.items()
    else:
        num_iters1, num_iters2 = iter(values1), iter(values2)

    min_, max_ = None, None
    for values in (num_iters1, num_iters2):
        for number in values:
            if min_ is None or min_ > number:
                min_ = number
            if max_ is None or max_ < number:
                max_ = number
    range_ = (min_, max_)

    #Now I calculate the distibution with the same range, in order to be able
    # to get the difference
    if low_memory:
        num_iters1, num_iters2 = values1.items(), values2.items()
    else:
        num_iters1, num_iters2 = iter(values1), iter(values2)

    distrib1 = _create_distribution(num_iters1,  range_=range_,
                                    low_memory=low_memory)
    distrib2 = _create_distribution(num_iters2,  range_=range_,
                                    low_memory=low_memory)

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

def general_seq_statistics(sequences, distrib_fhand=None, plot_fhand=None):
    '''Given a sequence iterator it calculates some general statistics.

    The statistics will be written into the given distrib_fhand (if given) and
    a dict with them will be returned.
    It calculates the total sequence length, the average sequence length, the
    total masked sequence length, and the number of sequences.
    '''
    stats   = {}
    stats['seq_length']      = 0
    stats['num_sequences']   = 0
    stats['masked_length']   = 0
    stats['max_length']      = None
    stats['min_length']      = None
    stats['length_variance'] = None
    seq_quality_average      = []

    seq_len_list = []
    first        = True
    for seq in sequences:
        seq_len = len(seq)
        stats['num_sequences'] += 1
        stats['seq_length']    += seq_len

        # variance, to calculate the variance of the sequences I need a list of
        # the sequences. I wild use it o calculate the max and min
        seq_len_list.append(seq_len)

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
                if qual is not None:
                    has_qual = True
                else:
                    has_qual = False
            except AttributeError:
                has_qual = False

        if has_qual:
            qual = seq.qual
            qual_average = sum(qual)/float(seq_len)
            seq_quality_average.append((seq_len, qual_average))

    # variance, max and min, average
    stats['length_variance']    = numpy.var(seq_len_list)
    stats['max_length']         = max(seq_len_list)
    stats['min_length']         = min(seq_len_list)
    stats['seq_length_average'] = numpy.average(seq_len_list)

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

def histogram(numbers, bins, range_=None):
    '''An alternative implemetation to the numpy.histogram.

    The main difference is that this implementation can use a CachedFileList to
    save memory
    '''
    #pylint: disable-msg=W0622
    #is this a FileCachedList?
    cached = False
    if 'items' in dir(numbers):
        cached = True
    #an iterator for the numbers
    if cached:
        num_iter = numbers.items()
    else:
        num_iter = iter(numbers)
    if range_ is None:
        #i'm not using the min and max functions to save one pass through the
        #data
        min_, max_ = None, None
        for number in num_iter:
            if min_ is None or min_ > number:
                min_ = number
            if max_ is None or max_ < number:
                max_ = number
    else:
        min_, max_ = range_[0], range_[1]

    #now we can calculate the bin egdes
    distrib_span = max_ - min_
    bin_span     = distrib_span / float(bins)
    bin_edges = [min_ + bin * bin_span for bin in range(bins + 1)]

    #now we calculate the distribution
    distrib = [0] * bins
    #an iterator for the numbers
    if cached:
        num_iter = numbers.items()
    else:
        num_iter = iter(numbers)
    for number in num_iter:
        if number == max_:
            #the last value go into the last bin
            bin = bins - 1
        else:
            bin = int(float(number - min_) / bin_span)
        distrib[bin] += 1

    return (distrib, bin_edges)
