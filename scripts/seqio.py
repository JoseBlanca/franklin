#!/usr/bin/env python
'It inputs and ouputs sequences, usually to change their format'

from optparse import OptionParser
from biolib.utils.seqio_utils import seqio

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-s', '--inseqfile', dest='inseqfile',
                    help='input sequence file')
    parser.add_option('-q', '--inqualfile', dest='inqualfile', default=None,
                      help='input quality file')
    parser.add_option('-f', '--informat', dest="informat", default=None,
                      help='input file format')
    parser.add_option('-t', '--outseqfile', dest='outseqfile',
                    help='output sequence file')
    parser.add_option('-r', '--outqualfile', dest='outqualfile', default=None,
                      help='output quality file')
    parser.add_option('-e', '--outformat', dest="outformat",
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
