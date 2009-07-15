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

from biolib.biolib_utils import draw_scatter

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

def _create_distribution(numbers, labels, distrib_fhand=None, plot_fhand=None):
    ''''Given a list of numbers it returns the distribution and it plots the
    histogram'''
    bins = 20
    #we do the distribution
    distrib, bin_edges = numpy.histogram(numbers, bins=bins, new=True)
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


def masked_seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None):
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
    labels = {'title': 'Masked sequence length distribution',
              'xlabel':'Masked length', 'ylabel':'Number of sequences'}
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand)

def seq_length_distrib(sequences, distrib_fhand=None, plot_fhand=None):
    '''It calculates the length distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    lengths = [len(seq) for seq in sequences]
    labels = {'title':'Sequence length distribution',
              'xlabel':'Sequence length', 'ylabel': 'Number of sequences'}
    return _create_distribution(lengths, labels, distrib_fhand, plot_fhand)

def seq_qual_distrib(sequences, distrib_fhand=None, plot_fhand=None):
    '''It calculates the qualities distribution for the given sequences.

    The results will be written in the given distrib_fhand.
    If plot_fhand is given a graphic with the distribution will be plotted.
    It returns the distribution and the bin_edges.
    '''
    #first we gather the lengths
    qualities = []
    for seq in sequences:
        for value in seq.qual:
            qualities.append(value)
    labels = {'title':'Sequence quality distribution', 'xlabel':'quality',
              'ylabel':'Number of bp'}
    return _create_distribution(qualities, labels, distrib_fhand, plot_fhand)
