#!/usr/bin/env python

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
It runs bwa.

It needs tha the database is indexed. It can run any type of sequence

Created on 08/10/2009

@author: peio
'''
from optparse import OptionParser
from biolib.biolib_seqio_utils import guess_seq_file_format
def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-s', '--seqs', dest='seqfile', help='Sequences file')
    parser.add_option('-f', '--format', dest='format',
                      help='Sequences file format')
    parser.add_option('-r', '--references', dest='reffiles',
                      help='reference seq file')
    parser.add_option('-i', '--refdb', dest='refdb', help='reference index')
    parser.add_option('-b', '--bamfile', dest='samfile', help='bam output file')

    return parser
def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.seqfile is None:
        parser.error('seqfile is mandatory')
    else:
        seqfhand = open(options.seqfile)

    if options.reffile is None:
        parser.error('file is mandatory')
    else:
        reffhand = open(options.seqfile)

    if options.format is None:
        file_format = guess_seq_file_format(seqfhand)
    else:
        file_format = options.format

    refdb = options.refdb
    if options.bamfile 
    


if __name__ == '__main__':
    main()