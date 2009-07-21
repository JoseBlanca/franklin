#!/usr/bin/env python

'Given a qual file it generates an histogram with the qualities distribution.'

from optparse import OptionParser

from Bio import SeqIO
import matplotlib.pyplot as plt

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser(version='0.1')
    parser.add_option('-i', '--infile', dest='infile', help='sequence file')
    return parser.parse_args()

def qualities(seqs):
    'Given a seq record iterator it returns a list with all qualities'
    quals = []
    for seq in seqs:
        quals.extend(seq.letter_annotations["phred_quality"])
    return quals

def main():
    'the main sub'
    options = parse_options()[0]
    fhand   = open(options.infile, 'r')
    seqs    = SeqIO.parse(fhand, 'qual')
    quals   = qualities(seqs)

    #ploting the figure
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.hist(quals, bins=20)
    axes.set_xlabel('quality')
    axes.set_ylabel('Num of bases')
    plt.show()


if __name__ == '__main__':
    main()
