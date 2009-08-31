#!/usr/bin/env python
'''It calculates a distribution for some sequence property.

It requires a sequence file (optionally with qualities). It calculates the
distribution of some property in the sequence files (like the sequence length).

The distribution types available are:
    - sequence length:        seq_length_distrib
    - quality:                qual_distrib
    - masked sequence length: masked_seq_distrib

The script can generate a plot with the distribution and/or a text file with the
distribution values.
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
from biolib.statistics import seq_distrib

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-s', '--seqfile', dest='seqfile1',
                      help='Seq1 input fasta')
    parser.add_option('-q', '--qualfile', dest='qualfile1',
                      help='Quality1 fasta file')
    parser.add_option('-f', '--seq_format', dest='format1',
                      help='sequence file format')
    parser.add_option('-p', '--plot', dest='plot',
                      help='plot output file')
    parser.add_option('-w', '--distrib', dest='distrib',
                      help='distribution text output file')
    parser.add_option('-t', '--type', dest='kind',
                      help='type of distribution')

    return parser.parse_args()

#pylint:disable-msg=R0912
def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''

    options = parse_options()[0]

    io_fhands = {}
    if options.kind is None:
        raise RuntimeError('Distribution type is mandatory (-t distrib_type)')
    else:
        kind = options.kind

    if options.seqfile1 is None:
        raise RuntimeError('Need seq file1')
    else:
        io_fhands['seqfile1'] = open(options.seqfile1, 'r')

    if options.qualfile1 is None:
        io_fhands['qualfile1'] = None
    else:
        io_fhands['qualfile1'] = open(options.qualfile1, 'r')

    if options.plot is None:
        io_fhands['plot'] = None
    else:
        io_fhands['plot'] = open(options.plot, 'w')

    if options.distrib is None:
        io_fhands['distrib'] = None
    else:
        io_fhands['distrib'] = open(options.distrib, 'w')

    return kind, io_fhands, options.format1

def main():
    'The main function'
    kind, io_fhands, format = set_parameters()

    if io_fhands['qualfile1']:
        seqs1 = seqs_in_file(io_fhands['seqfile1'], io_fhands['qualfile1'],
                             format=format)
    else:
        seqs1 = seqs_in_file(io_fhands['seqfile1'], format=format)

    seq_distrib(sequences=seqs1, kind=kind, distrib_fhand=io_fhands['distrib'],
                plot_fhand=io_fhands['plot'], low_memory=True)

if __name__ == '__main__':
    main()
