#!/usr/bin/env python

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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

'It inputs and ouputs sequences, usually to change their format'

from optparse import OptionParser
from franklin.utils.seqio_utils import seqio

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-s', '--inseqfile', dest='inseqfile',
                    help='input sequence file')
    parser.add_option('-q', '--inqualfile', dest='inqualfile', default=None,
                      help='input quality file')
    parser.add_option('-f', '--informat', dest="informat", default=None,
                      help='input file format')
    parser.add_option('-o', '--outseqfile', dest='outseqfile',
                    help='output sequence file')
    parser.add_option('-l', '--outqualfile', dest='outqualfile', default=None,
                      help='output quality file')
    parser.add_option('-t', '--outformat', dest="outformat",
                      help='output file format')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    infhand_seq = open(options.inseqfile)
    if options.inqualfile:
        infhand_qual = open(options.inqualfile)
    else:
        infhand_qual = None

    format_in = options.informat

    outfhand_seq = open(options.outseqfile, 'w')
    if options.outqualfile:
        outfhand_qual = open(options.outqualfile, 'w')
    else:
        outfhand_qual = None
    format_out = options.outformat

    return (infhand_seq, infhand_qual, format_in, outfhand_seq, outfhand_qual,
            format_out)


def main():
    'The main function'
    (in_seq_fhand, in_qual_fhand, in_format, out_seq_fhand, out_qual_fhand,
     out_format) = set_parameters()

    seqio(in_seq_fhand=in_seq_fhand, in_qual_fhand=in_qual_fhand,
          in_format=in_format,
          out_seq_fhand=out_seq_fhand, out_qual_fhand=out_qual_fhand,
          out_format=out_format)

if __name__ == '__main__':
    main()
