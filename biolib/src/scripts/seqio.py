#!/usr/bin/env python
'It inputs and ouputs sequences, usually to change their format'

from optparse import OptionParser
from biolib.biolib_seqio_utils import seqs_in_file, write_seqs_in_file
from Bio import SeqIO

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
    (infhand_seq, infhand_qual, format_in, outfhand_seq, outfhand_qual,
     format_out) = set_parameters()

    if  infhand_qual is not None or outfhand_qual is not None:
        seqs = seqs_in_file(seq_fhand=infhand_seq, qual_fhand=infhand_qual,
                        format=format_in, create_seqrecord=True)
        write_seqs_in_file(seqs, seq_fhand=outfhand_seq,
                           qual_fhand=outfhand_qual,format=format_out)
    else:
        SeqIO.convert(infhand_seq, format_in, outfhand_seq, format_out)


if __name__ == '__main__':
    main()