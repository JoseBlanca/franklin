#!/usr/bin/env python
'''

Created on 2009 eka 2

@author: peio
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
from biolib.biolib_utils import draw_scatter
from biolib.alignment_search_result import (BlastParser,
                                            alignment_results_scores)
import numpy
import pylab

def main():
    '''The main section'''
    parser = OptionParser('usage: %prog -i blast.xml', version='%prog 0.1')
    parser.add_option('-i', '--infile', dest='infile', help='blast xml')
    msg = 'sum hits in the distribution, not lengths'
    parser.add_option('-t', '--hits', dest='use_length',
                      action='store_true', default=False, help=msg)
    parser.add_option('-c', '--incompat', dest='do_incompat',
                      action='store_true', default=False, help=msg)
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script at least needs an input file (blast xml output)')
    else:
        infile = options.infile
    #use_length = options.use_length

    bins = 20

    blasts = BlastParser(fhand=open(infile, 'r'))
    #the values for the distribution
    score_keys = ['similarity']
    if options.do_incompat:
        score_keys.append('d_incompatibility')
    scores = alignment_results_scores(blasts, score_keys)
    #the distribution
    if options.do_incompat:
        #distrib, x_edges, y_edges = numpy.histogram2d(scores[0], scores[1],
        #                                              bins=bins)
        distrib = numpy.histogram2d(scores[0], scores[1], bins=bins)[0]
    else:
        distrib, bin_edges = numpy.histogram(scores, bins=bins)
    #the drawing
    if options.do_incompat:
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
        draw_scatter(x_axe=bin_edges[:-1], y_axe=distrib)
    return

if __name__ == '__main__':
    main()
