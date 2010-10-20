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

'''
it removes annotation from seqs in the file
'''

from optparse import OptionParser
import sys
from franklin.seq.readers import guess_seq_file_format, seqs_in_file
from franklin.seq.writers import write_seqs_in_file

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--inseqfile', dest='infile',
                    help='input sequence file')
    parser.add_option('-o', '--outseqfile', dest='outfile',
                    help='output sequence file')
    parser.add_option('-r', '--remove_annotations', dest='rm_annot',
                      help='Annotation to remove. comma separated')

    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error = 'Infile needed'
    else:
        infhand = open(options.infile)

    if options.outfile is None:
        outfhand = sys.stdout
    else:
        outfhand = open(options.outfile, 'w')

    if options.rm_annot is None:
        parser.error = 'Need annotation name to remove'
    else:
        rm_annots = options.rm_annot.split(',')

    return infhand,  outfhand, rm_annots

def main():
    'The main'
    # get parameters
    infhand, outfhand, rm_annots = set_parameters()

    # guess file format
    format_ = guess_seq_file_format(infhand)

    #remove annotations
    seqs = remove_annotation(infhand, format_, rm_annots)

    # write seqs in file
    write_seqs_in_file(seqs, seq_fhand=outfhand, format=format_)

def remove_annotation(in_fhand, format_, rm_annots):
    'it remoes annotations and yields the seqs'
    for seq in seqs_in_file(in_fhand, format=format_):
        for rm_annot in rm_annots:
            seq.remove_annotations(rm_annot)
        yield seq

if __name__ == '__main__':
    main()
