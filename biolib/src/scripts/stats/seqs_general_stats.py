#!/usr/bin/env python
'''This script calculates some general statistics for a sequence file.

It calculates:
    - Total sequence length
    - Average sequence length
    - Masked sequence length
    - Number of sequences
    - Quality Average
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
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

from optparse import OptionParser
from biolib.biolib_seqio_utils import seqs_in_file
from biolib.statistics import general_seq_statistics


def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-s', '--fastafile', dest='fastafile',
                      help='Sanger input fasta')
    parser.add_option('-q', '--qualfile', dest='qualfile',
                      help='Quality fasta file')
    parser.add_option('-w', '--resultfile', dest='resultfile',
                      help='result output file')

    return parser.parse_args()

def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''
    options = parse_options()[0]
    if options.fastafile is None:
        raise RuntimeError('Script at least needs an input file (fasta)')
    else:
        fhand_seq = open(options.fastafile, 'r')

    if options.qualfile is None:
        fhand_qual = None
    else:
        fhand_qual = options.qualfile

    if options.resultfile is None:
        result_file = None
    else:
        result_file = open(options.resultfile, 'w')

    return fhand_seq, fhand_qual, result_file


def main():
    'Main section'
    fhand_seq, fhand_qual, result_file = set_parameters()
    seqs = seqs_in_file(fhand_seq, fhand_qual)

    stats = general_seq_statistics(seqs, distrib_fhand=result_file)

    for key, value in stats.items():
        print '%-19s : %d' % (key, value)


if __name__ == '__main__':
    main()
