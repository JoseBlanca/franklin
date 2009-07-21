#!/usr/bin/env python
'''
This script calculates the distribution difference between two seqs applying to
it the different distribution calculation method that we have. It plots the
differene
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
from biolib.biolib_utils import seqs_in_file
from biolib.statistics import seq_distrib_diff

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-s', '--seqfile1', dest='seqfile1',
                      help='Seq1 input fasta')
    parser.add_option('-q', '--qualfile1', dest='qualfile1',
                      help='Quality1 fasta file')
    parser.add_option('-l', '--seqfile2', dest='seqfile2',
                      help='Seq2 input fasta')
    parser.add_option('-k', '--qualfile2', dest='qualfile2',
                      help='Quality2 fasta file')
    parser.add_option('-p', '--plot', dest='plot',
                      help='plot output file')
    parser.add_option('-w', '--distrib', dest='distrib',
                      help='distribution text output file')
    parser.add_option('-t', '--type', dest='kind',
                      help='type of analisis')

    return parser.parse_args()

#pylint:disable-msg=R0912
def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''

    options = parse_options()[0]

    io_fhands = {}
    if options.kind is None:
        raise RuntimeError('Stat analisis type is mandatory, use -t stat_type')
    else:
        kind = options.kind

    if options.seqfile1 is None:
        raise RuntimeError('Need seq file1')
    else:
        io_fhands['seqfile1'] = open(options.seqfile1, 'r')

    if options.qualfile1 is None:
        if kind in ['qual_distrib']:
            raise RuntimeError('Need first qual file to calculate %s analisis'\
                                % kind)
        else:
            io_fhands['qualfile1'] = None
    else:
        io_fhands['qualfile1'] = open(options.qualfile1, 'r')

    if options.seqfile2 is None:
        raise RuntimeError('Need second seq file')
    else:
        io_fhands['seqfile2'] = open(options.seqfile2, 'r')

    if options.qualfile2 is None:
        if kind in ['qual_distrib']:
            raise RuntimeError('Need second qual file to calculate %s analisis'\
                                % kind)
        else:
            io_fhands['qualfile2'] = None
    else:
        io_fhands['qualfile2'] = open(options.qualfile2, 'r')

    if options.plot is None:
        io_fhands['plot'] = None
    else:
        io_fhands['plot'] = open(options.plot, 'w')

    if options.distrib is None:
        io_fhands['distrib'] = None
    else:
        io_fhands['distrib'] = open(options.distrib, 'w')

    return kind, io_fhands

def main():
    'The main function'
    kind, io_fhands = set_parameters()

    seqs1 = seqs_in_file( io_fhands['seqfile1'],  io_fhands['qualfile1'])
    seqs2 = seqs_in_file( io_fhands['seqfile2'],  io_fhands['qualfile2'])

    seq_distrib_diff(seqs1, seqs2, kind, distrib_fhand=io_fhands['distrib'],
                     plot_fhand=io_fhands['plot'])

if __name__ == '__main__':
    main()
