'It inputs and ouputs sequences, usually to change their format'

from optparse import OptionParser
import Bio
from biolib.biolib_seqio_utils import seqs_in_file, write_seqs_in_file

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
    return parser.parse_args()

def main():
    'The main function'
    # Set parameters
    opts = parse_options()[0]
    in_seq_fhand = open(opts.inseqfile)
    if opts.inqualfile:
        in_qual_fhand = open(opts.inqualfile)
    else:
        in_qual_fhand = None
    in_format = opts.informat

    out_seq_fhand = open(opts.outseqfile, 'w')
    if opts.outqualfile:
        out_qual_fhand = open(opts.outqualfile, 'w')
    else:
        out_qual_fhand = None
    out_format = opts.outformat

    seqs = seqs_in_file(seq_fhand=in_seq_fhand,
                        qual_fhand=in_qual_fhand,
                        format=in_format,
                        create_seqrecord=True)
    write_seqs_in_file(seqs, seq_fhand=out_seq_fhand,
                       qual_fhand=out_qual_fhand,
                       format=out_format)

if __name__ == '__main__':
    main()