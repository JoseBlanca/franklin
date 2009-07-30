#!/usr/bin/env python
'''
This script takes sequences from a fasta file. Takes each sequence and removes
the masked section. The resulting sequence sections are returned as new
sequences

Created on 2009 uzt 24

@author: peio
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
from biolib.seq_cleaner import  split_seq_by_masked_regions
from biolib.biolib_seqio_utils import seqs_in_file, write_fasta_file

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-s', '--fastafile', dest='fastafile',
                      help='Sanger input fasta')
    parser.add_option('-q', '--qualfile', dest='qualfile',
                      help='Quality fasta file')
    parser.add_option('-o', '--seqoutput', dest='out_seq',
                      help='sequence output file')
    parser.add_option('-u', '--qualoutput', dest='out_qual',
                      help='quality output file')
    parser.add_option('-m', '--min_length', dest='minlength',
                      type="int", help='Minimun sequence length for new seqs')
    return parser.parse_args()

def set_parameters():
    '''It sets the parameters for the script.

    For some options there are some default values that will be applied here.
    '''
    options = parse_options()[0]

    ############### input output parameters ################
    io_fhands = {}

    if options.fastafile is None:
        raise RuntimeError('Script at least needs an input file (fasta)')
    else:
        io_fhands['in_seq'] = open(options.fastafile, 'r')

    if options.qualfile is None:
        io_fhands['in_qual'] = None
    else:
        io_fhands['in_qual'] = open(options.qualfile, 'r')

    if options.out_seq is None:
        io_fhands['out_seq'] = open('result.seq.fasta', 'w')
    else:
        io_fhands['out_seq'] = open(options.out_seq, 'w')

    if options.out_qual is None:
        if  io_fhands['in_qual'] is None:
            io_fhands['out_qual'] = None
        else:
            io_fhands['out_qual'] = open('result.qual.fasta', 'w')
    else:
        io_fhands['out_qual'] = open(options.out_qual, 'w')
    ##########################################################
    if options.minlength is not None:
        minlength = options.minlength
    else:
        raise RuntimeError('Minimun resulting length seq needed')

    return io_fhands,  minlength



def main():
    'The main part of the script'
    io_fhands, minlength = set_parameters()

    #Get sequences from input files
    seq_iter     = seqs_in_file(io_fhands['in_seq'], io_fhands['in_qual'])

    # split new long seqs

    new_seq_iter = split_seq_by_masked_regions(seq_iter, minlength)

    # Write cutted seqs to a new fasta
    write_fasta_file(new_seq_iter, io_fhands['out_seq'], io_fhands['out_qual'])

if __name__ == '__main__':
    #main()
    import cProfile
    cProfile.run('main()')

