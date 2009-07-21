#!/usr/bin/env python

'''Given a sequence file it generates an histogram with the sequence length
distribution.'''

from optparse import OptionParser

from Bio import SeqIO
import matplotlib.pyplot as plt

from biolib.biolib_utils import guess_seq_file_format

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser(version='0.1')
    parser.add_option('-i', '--infile', dest='infile', help='sequence file')
    return parser.parse_args()

def main():
    'the main sub'
    options = parse_options()[0]

    fhand   = open(options.infile, 'r')
    format_ = guess_seq_file_format(fhand)
    seqs    = SeqIO.parse(fhand, format_)
    lengths = [len(seq) for seq in seqs]

    #ploting the figure
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.hist(lengths, bins=20)
    axes.set_xlabel('length')
    axes.set_ylabel('Num of sequences')
    plt.show()


if __name__ == '__main__':
    main()
