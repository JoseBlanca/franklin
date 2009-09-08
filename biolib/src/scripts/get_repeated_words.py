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
# along with biolib. If not, see <http://www.gnu.org/licenses/>.
'''
Created on 04/09/2009

@author: peio
'''
from optparse import OptionParser
from biolib.biolib_seqio_utils import seqs_in_file
from biolib.statistics import create_distribution
from biolib.words import cluster_seqs_by_words, filter_low_abundant_words

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -s seqfile, -w wordsize, -o outputfile')
    #I/O options
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-s', '--infile', dest='infile', help='input file')
    parser.add_option('-f', '--format', dest='format', help='file format')

    #Script parameters
    parser.add_option('-l', '--length', dest='wsize', help='Word length',
                      type='int', default=13)
    parser.add_option('-w', '--word_abundance', dest='word_abundance',
                      help='filter low abundance words', type='int', default=50)

    #stats options
    parser.add_option('-p', '--plotfile', dest='plotfile', help='Output file')
    parser.add_option('-d', '--distrib', dest='distrib',
                      help='create an abundance distribution')

    return parser

def set_parameters():
    '''It sets the parameters for the script.

    For some options there are some default values that will be applied here.
    '''
    parser = parse_options()
    options = parser.parse_args()[0]
    #I/O options
    if options.infile is None:
        parser.error('Script at least needs an input file')
    else:
        in_fhand = open(options.infile)

    if options.outfile is None:
        out_fhand = open(options.infile + '.out', 'w')
    else:
        out_fhand = open(options.outfile, 'w')
    if options.format is None:
        format_ = 'fasta'
    else:
        format_ = options.format
    #Script parameters
    wsize          = options.wsize
    word_abundance = options.word_abundance

    #stats options
    if options.plotfile is None:
        plot_fhand = None
    else:
        plot_fhand = open(options.plotfile, 'w')

    if options.distrib is None:
        distrib_fhand = None
    else:
        distrib_fhand = open(options.distrib, 'w')

    return (in_fhand, out_fhand, format_, wsize, word_abundance, distrib_fhand,
            plot_fhand)

def main():
    'The main function'
    #Set scripts input parameters
    (in_fhand, out_fhand, format_, wsize, word_abundance, distrib_fhand,
     plot_fhand) = set_parameters()

    #Get seqs from file
    seqs = seqs_in_file(in_fhand, format=format_)

    # Count words in seqs
    cluster_seqs = cluster_seqs_by_words(seqs, wsize, True)[0]

    # Stats: only if it is asked
    if distrib_fhand is not None or plot_fhand is not None:
        #we cache the result in disk
        create_distribution(cluster_seqs.non_empty_values(), plot_fhand=plot_fhand,
                            distrib_fhand=distrib_fhand, range_= (1, None),
                            labels={'title':'word distribution',
                                    'xlabel':'times found',
                                    'ylabel':'number_of_words'},
                            low_memory=True)

    #Filter the less frecuent words. We get most repeated words
    words = filter_low_abundant_words(cluster_seqs, word_abundance, wsize)
    for word_abundance in words:
        word, abundance = word_abundance
        out_fhand.write('%s: %i\n' % (word, abundance))

        
if __name__ == '__main__':
    main()
