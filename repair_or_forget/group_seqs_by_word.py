# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

'''
This scripts groups the sequences looking at the words they have.

@author: peio
'''
from optparse import OptionParser
import sys, itertools
from franklin.biolib_seqio_utils import seqs_in_file
from franklin.statistics import create_distribution
from franklin.words import cluster_seqs_by_words, filter_low_abundant_words

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -s seqfile, -l wordsize, -o outputfile')
    #I/O options
    parser.add_option('-o', '--outfile', dest='outfile', help='Output file')
    parser.add_option('-s', '--infile', dest='infile', help='input file')
    parser.add_option('-f', '--format', dest='format', help='file format')

    #Script parameters
    parser.add_option('-l', '--length', dest='wsize', help='Word length',
                      type='int', default=13)

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

    return (in_fhand, out_fhand, format_, wsize)

def write_cluster_file(cluster_seqs, seq_names, out_fhand):
    '''This function writes the cluster info to a file to be used by a
    clustering program'''
    for cluster in cluster_seqs.non_empty_values():
        if cluster is None or len(cluster) < 2:
            continue
        first = seq_names[cluster.pop()]
        for i in range(len(cluster)):
            out_fhand.write("%s %s\n" % (first, seq_names[cluster.pop()]))

def main():
    'The main function'
    #Set scripts input parameters
    (in_fhand, out_fhand, format_, wsize) = set_parameters()

    #Get seqs from file
    seqs = seqs_in_file(in_fhand, format=format_)

    # Group seqs by words
    cluster_seqs, seq_names = cluster_seqs_by_words(seqs, wsize)

    #write clusters in pairs to a file
    write_cluster_file(cluster_seqs, seq_names, out_fhand)

if __name__ == '__main__':
    main()
