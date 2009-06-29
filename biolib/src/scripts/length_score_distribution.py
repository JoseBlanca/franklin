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
                                            generate_score_distribution)

def main():
    '''The main section'''
    parser = OptionParser('usage: %prog -i blast.xml', version='%prog 0.1')
    parser.add_option('-i', '--infile', dest='infile', help='blast xml')
    msg = 'sum lengths in the distribution, not hits'
    parser.add_option('-l', '--length', dest='use_length',
                      action='store_false', default=True, help=msg)
    parser.add_option('-c', '--incompat', dest='do_incompat',
                      action='store_true', default=False, help=msg)
    (options, args) = parser.parse_args()

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    use_length = options.use_length
 
    blasts = BlastParser(fhand=open(infile, 'r'))
    distrib = generate_score_distribution(results=blasts,
                                          score_key  = 'similarity',
                                          nbins      = 20,
                                          use_length = use_length,
                                calc_incompatibility = options.do_incompat)
    print('distribution -> ' + str(distrib))
    if options.do_incompat:
        #draw 3d distrib
        dis   = distrib['distribution']
        x_axe = distrib['similarity_bins']
        y_axe = distrib['incompatibility_bins']
        z_values = []
        for x_index, score_list in enumerate(dis):
            x_value = x_axe[x_index]
            x_values = []
            for y_index, z_value in enumerate(score_list):
                y_value = y_axe[y_index]
                if z_value is None:
                    z_value = 0
                x_values.append(z_value)
            z_values.append(x_values)
        #the drawing
        import pylab
        axes = pylab.subplot(111)
        image = pylab.imshow(z_values)
        image.set_interpolation('bilinear')
        pylab.show()
    else:
        draw_scatter(x_axe=distrib['bins'][:-1], y_axe=distrib['distribution'])

if __name__ == '__main__':
    main()
