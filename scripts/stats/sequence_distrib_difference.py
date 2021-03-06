#!/usr/bin/env python
'''It calculates the difference between two distributions.

It requires to sequence files (optionally with qualities). It calculates the
distribution of some property in the sequence files (like the sequence length)
and it creates a new distribution subtracting the one obtained for the first
sequence file with the one obtained with the second sequence file.

The distribution types available are:
    - sequence length:        seq_length_distrib
    - quality:                qual_distrib
    - masked sequence length: masked_seq_distrib

The script can generate a plot with the distribution and/or a text file with the
distribution values.
'''

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

from optparse import OptionParser
from franklin.utils.seqio_utils import seqs_in_file
from franklin.statistics import seq_distrib_diff

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
    help1  = 'type of distribution:masked_seq_distrib,seq_length_distrib,'
    help1 += 'qual_distrib'

    parser.add_option('-t', '--type', dest='kind', help=help1)

    parser.add_option('-f', '--inputformat1', dest="file_format1",
     help='Input file format: fasta(default), fastq, fastaq, fastq-solexa, ...')

    parser.add_option('-g', '--inputformat2', dest="file_format2",
     help='Input file format: fasta(default), fastq, fastaq, fastq-solexa, ...')

    return parser.parse_args()

#pylint:disable-msg=R0912
def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''

    options = parse_options()[0]
    ########### file format ###############################
    if options.file_format1 is None:
        input_file_format1 = None
    else:
        input_file_format1 = options.file_format1

    if options.file_format2 is None:
        input_file_format2 = None
    else:
        input_file_format2 = options.file_format2
    ##########################################################


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

    if options.seqfile2 is None:
        raise RuntimeError('Need second seq file')
    else:
        io_fhands['seqfile2'] = open(options.seqfile2, 'r')

    if options.qualfile2 is None:
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



    return kind, io_fhands, input_file_format1, input_file_format2

def main():
    'The main function'
    kind, io_fhands, input_file_format1, input_file_format2 = set_parameters()

    seqs1 = seqs_in_file( io_fhands['seqfile1'],  io_fhands['qualfile1'],
                          input_file_format1)
    seqs2 = seqs_in_file( io_fhands['seqfile2'],  io_fhands['qualfile2'],
                          input_file_format2)

    seq_distrib_diff(seqs1, seqs2, kind, distrib_fhand=io_fhands['distrib'],
                     plot_fhand=io_fhands['plot'])

if __name__ == '__main__':
    main()
