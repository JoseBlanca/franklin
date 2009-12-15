#!/usr/bin/env python
'''It creates scores distributions for the blast files.

Given a blast file this script is capable of creating two types of scoring
distributions. In the most simple one we get a plot of the number of blast
hits against the similarity. For every similarity percentage we get how many
blast hits are in the file.

For the hits besides the similarity we can calculate how many bases that should
be aligned between the query and the subject are not aligned. E.g.
query    ---------------->
         |||||||
subject  ---------------->
                <--------> Not aligned region
We can plot the distribution of hits taking into account the similarity
percentage and the percentage of the not aligned region. In that case we would
get a 3-D distribution plot. To get that you should chose the --incompat option.
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with Foobar. If not, see <http://www.gnu.org/licenses/>.

from optparse import OptionParser
from biolib.statistics import draw_scatter
from biolib.alignment_search_result import (BlastParser,
                                                       alignment_results_scores)
from biolib.statistics import create_distribution

import numpy
import pylab

def parse_options():
    'It parses the command line arguments'
    '''The main section'''
    parser = OptionParser('usage: %prog -i blast.xml', version='%prog 0.1')
    parser.add_option('-i', '--infile', dest='infile', help='blast xml')
    parser.add_option('-o', '--outfile', dest='outfile', help='png output')
    msg = 'sum hits in the distribution, not lengths'
    parser.add_option('-t', '--hits', dest='use_length',
                      action='store_true', default=False, help=msg)
    parser.add_option('-c', '--incompat', dest='do_incompat',
                      action='store_true', default=False, help=msg)
    parser.add_option('-l', '--low_memory', dest='low_memory',
                      action='store_true', default=False, help='Use low_memory')

    return parser

def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''
    parser  = parse_options()
    options = parser.parse_args()[0]
    if options.infile is None:
        parser.error('Script at least needs an input file (blast xml output)')
    else:
        infhand = open(options.infile)
    if options.outfile is None:
        outfhand = None
    else:
        outfhand = open(options.outfile, 'w')

    return infhand, outfhand, options.do_incompat, options.low_memory

def main():
    '''The main section'''

    # Get parameters
    infhand, outfhand, do_incompat, low_memory = set_parameters()
    bins = 20
    range_ = (95, 100)

    # Parse blast results
    blasts = BlastParser(infhand)

    # The values for the distribution
    score_keys = ['similarity']

    if do_incompat:
        score_keys.append('d_incompatibility')
    scores = alignment_results_scores(blasts, score_keys)

    # The distribution
    if do_incompat:
        #distrib, x_edges, y_edges = numpy.histogram2d(scores[0], scores[1],
        #                                              bins=bins)
        distrib = numpy.histogram2d(scores[0], scores[1], bins=bins)[0]
    else:
        result = create_distribution(scores, range=range_, bins=bins,
                                     low_memory=low_memory)
        distrib   = result['distrib'][0]
        bin_edges = result['bin_edges'][0]

    # The drawing
    if do_incompat:
        #fig = pylab.figure()
        pylab.figure()
        #axes = Axes3D(fig)
        #axes.plot_surface(x_edges[:-1], y_edges[:-1], distrib)
        #axes = pylab.subplot(111)
        pylab.subplot(111)
        image = pylab.imshow(distrib)
        image.set_interpolation('bilinear')
        pylab.show()
    else:
        draw_scatter(x_axe=bin_edges[:-1], y_axe=distrib, fhand=outfhand)
    return

if __name__ == '__main__':
    main()
