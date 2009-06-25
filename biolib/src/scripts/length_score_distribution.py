# -*- coding= UTF-8
'''
Created on 2009 eka 2

@author: peio
'''
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
    (options, args) = parser.parse_args()

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
    use_length = options.use_length
 
    blasts = BlastParser(fhand=open(infile, 'r'))
    distrib = generate_score_distribution(results=blasts,
                                          score_key='similarity',
                                          nbins    = 20,
                                          use_length   = use_length)
    print('distribution -> ' + str(distrib))
    draw_scatter(x_axe=distrib['bins'][:-1], y_axe=distrib['distribution'])

if __name__ == '__main__':
    main()
